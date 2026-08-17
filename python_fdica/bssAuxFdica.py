"""補助関数法による周波数領域独立成分分析（AuxFDICA）の本体。

複数マイクで観測した混合音をSTFTで周波数領域へ移し、周波数ごとの
分離行列を反復更新して音源を分離する。分離後はprojection backで振幅を
基準マイクへ合わせ、置換問題を解いてから時間波形へ戻す。
MATLAB版 ``bssAuxFdica.m`` と対応する関数構成を保っている。

このファイルで頻出するテンソル軸は次の通りである。
``sample`` は時間サンプル、``frequency`` は周波数ビン、``frame`` は
STFTの時間フレーム、``channel`` はマイク、``source`` は推定音源を表す。
"""

from __future__ import annotations

import time
from collections.abc import Sequence

import torch
import torch.nn.functional as torch_functional

from .permSolverCor import permSolverCor
from .permSolverDoa import permSolverDoa
from .permSolverIps import permSolverIps
from .stft import dgt_istft, dgt_stft


def local_whitening(X: torch.Tensor, nSrc: int) -> torch.Tensor:
    """周波数ごとにPCA白色化し、成分間の相関とスケール差を取り除く。

    入力 ``X`` は ``(frequency, frame, channel)``、出力 ``Y`` は
    ``(frequency, frame, source)``。白色化後の共分散行列は単位行列に
    近くなり、FDICAが分離行列を求めやすくなる。
    """
    nFreq, nFrame, _ = X.shape
    Y = torch.empty((nFreq, nFrame, nSrc), dtype=X.dtype, device=X.device)
    eps = torch.finfo(X.real.dtype).eps
    for iFreq in range(nFreq):
        # 1周波数分を、線形代数で一般的な (channel, frame) へ転置する。
        Xi = X[iFreq].mT
        # Xi Xi^H / フレーム数は、チャンネル間の複素共分散行列である。
        covariance = Xi @ Xi.mH / nFrame
        # Hermite行列専用のeighを使い、固有値と直交固有ベクトルを求める。
        eigenvalues, eigenvectors = torch.linalg.eigh(covariance)
        # 固有値が大きい順にnSrc成分を残す。これは主成分分析に相当する。
        idx = torch.argsort(eigenvalues, descending=True)[:nSrc]
        dP = eigenvectors[:, idx]
        # 主成分へ射影し、各成分を固有値の平方根で割って分散を1へそろえる。
        Y[iFreq] = ((dP.mH @ Xi) / torch.sqrt(eigenvalues[idx].clamp_min(eps))[:, None]).mT
    return Y


def local_calcFdicaCost(Y: torch.Tensor, W: torch.Tensor, srcModel: str) -> torch.Tensor:
    """FDICAの目的関数（コスト関数）を計算する。

    ``Y`` は ``(frequency, frame, source)``、``W`` は周波数ごとの分離行列。
    音源モデルに応じた罰則と、分離行列が退化しないためのlog行列式項から
    構成される。値そのものより、反復に伴う変化を見るために利用する。
    """
    nFrame = Y.shape[1]
    # log(0)を避けるため、行列式の絶対値を機械イプシロン以上に制限する。
    logdet = torch.log(torch.abs(torch.linalg.det(W)).clamp_min(torch.finfo(Y.real.dtype).eps)).sum()
    if srcModel == "LAP":
        return torch.abs(Y).sum() - 2 * nFrame * logdet
    return torch.log((torch.abs(Y) ** 2).clamp_min(torch.finfo(Y.real.dtype).eps)).sum() - 2 * nFrame * logdet


def local_plotSpectrogram(
    signal: torch.Tensor,
    sampFreq: float,
    fftSize: int,
    shiftSize: int,
    *,
    title: str,
    dynamicRange: float = 80.0,
    trunc: float = 0.0,
    normalize: bool = True,
) -> None:
    """MATLAB版DGTtoolに似たスペクトログラム・スペクトル・波形を描く。

    1チャンネル（または1音源）につき1枚のFigureを作る。左側に全区間の
    周波数スペクトル、右上に時間周波数スペクトログラム、右下に時間波形を
    配置する。描画は結果確認用なので、計算グラフから切り離してCPUへ移す。
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("isDraw=True requires matplotlib. Install it with `pip install matplotlib`.") from exc

    if signal.ndim != 2:
        raise ValueError("signal must be shaped (sample, channel/source)")
    # detachで勾配追跡を外し、matplotlibが扱えるCPUテンソルへ移動する。
    x = signal.detach().cpu()
    nSample, nSignal = x.shape
    window = torch.blackman_window(fftSize, periodic=True, dtype=x.dtype)
    # 窓関数による振幅低下を補正し、MATLAB版の表示スケールへ近づける。
    normConst = torch.sum(window).item() / 2 if normalize else 1.0
    spec = dgt_stft(x, fftSize, shiftSize).detach().cpu() / normConst
    xForFft = x / normConst
    eps = torch.finfo(spec.real.dtype).eps

    if sampFreq == 1:
        # サンプリング周波数1は、実時間ではなく正規化単位で表示する場合。
        fsPlot = 1.0
        freqUnit = "[periods/sample]"
        timeUnit = "[samples]"
        time = torch.arange(spec.shape[1], dtype=x.real.dtype) * shiftSize
        waveTime = torch.arange(nSample, dtype=x.real.dtype)
        xLimit = (0, max(nSample - 1, 1))
    else:
        # 通常の音声では、横軸を秒、縦軸をkHzへ変換して表示する。
        fsPlot = sampFreq / 1000
        freqUnit = "[kHz]"
        timeUnit = "[s]"
        time = torch.arange(spec.shape[1], dtype=x.real.dtype) * shiftSize / sampFreq
        waveTime = torch.arange(nSample, dtype=x.real.dtype) / sampFreq
        xLimit = (0, max((nSample - 1) / sampFreq, 1 / sampFreq))
    freq = torch.linspace(0, fsPlot / 2, spec.shape[0], dtype=x.real.dtype)

    for index in range(nSignal):
        # 振幅をdBへ変換する。epsでlog10(0)による-infを防ぐ。
        spectrum = 20 * torch.log10(torch.abs(torch.fft.rfft(xForFft[:, index])).clamp_min(eps))
        specDb = 20 * torch.log10(torch.abs(spec[:, :, index]).clamp_min(eps))
        vmax = float(specDb.max().item()) - trunc
        vmin = vmax - dynamicRange

        # MATLAB版と同様の領域配置をadd_gridspecで組み立てる。
        fig = plt.figure()
        fig.suptitle(f"{title} {index + 1}")
        grid = fig.add_gridspec(10, 14, wspace=0.0, hspace=0.0)
        axSpectrum = fig.add_subplot(grid[:8, :2])
        axSpec = fig.add_subplot(grid[:8, 2:])
        axWave = fig.add_subplot(grid[8:, 2:])
        axColor = fig.add_subplot(grid[8:, :1])

        specFreq = torch.linspace(0, fsPlot / 2, spectrum.numel(), dtype=x.real.dtype)
        axSpectrum.plot(spectrum.numpy(), specFreq.numpy())
        axSpectrum.set_xlim(float(spectrum.max().item()) - dynamicRange - trunc, float(spectrum.max().item()) - trunc)
        axSpectrum.set_ylim(0, fsPlot / 2)
        axSpectrum.invert_xaxis()
        axSpectrum.set_ylabel(f"Frequency {freqUnit}", fontsize=12)
        axSpectrum.tick_params(labelsize=10)

        image = axSpec.imshow(
            specDb.numpy(),
            origin="lower",
            aspect="auto",
            extent=(float(time[0]), float(time[-1]) if time.numel() else 0, 0, fsPlot / 2),
            vmin=vmin,
            vmax=vmax,
        )
        axSpec.set_xlim(*xLimit)
        axSpec.set_ylim(0, fsPlot / 2)
        axSpec.axis("off")

        axWave.plot(waveTime.numpy(), x[:, index].numpy())
        axWave.set_xlim(*xLimit)
        axWave.set_xlabel(f"Time {timeUnit}", fontsize=12)
        axWave.set_yticks([0])
        axWave.axhline(0, color="0.3", linewidth=0.8)
        axWave.tick_params(labelsize=10)

        fig.colorbar(image, cax=axColor)
        axColor.set_ylabel("Power [dB]", fontsize=11)
        axColor.tick_params(axis="x", bottom=False, labelbottom=False)


def local_plotCost(cost: torch.Tensor, nIter: int) -> None:
    """0回目から最終反復までのコスト関数の軌跡を折れ線で描く。"""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("isDraw=True requires matplotlib. Install it with `pip install matplotlib`.") from exc

    costCpu = cost.detach().cpu()
    fig, ax = plt.subplots()
    ax.plot(range(nIter + 1), costCpu[: nIter + 1].numpy())
    ax.set_xlabel("Number of iterations")
    ax.set_ylabel("Value of cost function")
    ax.grid(True)
    ax.tick_params(labelsize=12)


def local_showPlots() -> None:
    """作成した全Figureを表示し、利用者が閉じるまで処理を待つ。"""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("isDraw=True requires matplotlib. Install it with `pip install matplotlib`.") from exc

    plt.show(block=True)


def local_auxFdica(
    X: torch.Tensor,
    nIter: int,
    srcModel: str,
    isDraw: bool = False,
    verbose: bool = False,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """LAPまたはTVG音源モデルを使ってAuxFDICAの反復更新を実行する。

    入力 ``X`` は ``(frequency, frame, channel)``。各周波数にある分離行列
    ``W`` を反復的に更新し、推定スペクトログラム ``Y`` と分離行列、
    コスト履歴を返す。determined条件ではチャンネル数と音源数が等しい。
    """
    nFreq, nFrame, nCh = X.shape
    # 最初は「入力をそのまま出力する」単位行列を全周波数へ用意する。
    W = torch.eye(nCh, dtype=X.dtype, device=X.device).repeat(nFreq, 1, 1)
    Y = X.clone()
    # MATLAB版と同様、描画するときだけ初期値を含むコストを計算する。
    cost = torch.zeros(nIter + 1 if isDraw else nIter, dtype=X.real.dtype, device=X.device)
    if isDraw:
        cost[0] = local_calcFdicaCost(Y, W, srcModel)
    # 0除算や極端に大きな更新を避けるための、dtypeに応じた下限値。
    threshold = 10000 * torch.finfo(X.real.dtype).eps
    eye = torch.eye(nCh, dtype=X.dtype, device=X.device)
    for iIter in range(nIter):
        iterationStart = time.perf_counter()
        if verbose:
            print(f"[FDICA] 学習 {iIter + 1}/{nIter} 開始", flush=True)
        # LAPは振幅、TVGはパワーを補助変数の重みに使う。
        radius = torch.abs(Y) if srcModel == "LAP" else torch.abs(Y) ** 2
        invRadius = radius.clamp_min(threshold).reciprocal()
        for n in range(nCh):
            for f in range(nFreq):
                # 音源n・周波数fに対応する重み付き共分散行列Vnを作る。
                Xf = X[f].mT  # (channel, frame)
                Vn = (Xf * invRadius[f, :, n][None, :]) @ Xf.mH / nFrame
                # 無音に近い帯域や強く相関したマイク信号ではVnが特異行列に
                # なり得る。微小値を対角へ加え、数値的に解ける状態へする。
                meanPower = torch.real(torch.trace(Vn)) / nCh
                ridge = 100 * torch.finfo(X.real.dtype).eps * meanPower.clamp_min(1.0)
                VnReg = Vn + ridge * eye
                system = W[f] @ VnReg
                try:
                    # iterative projectionの更新式 (W Vn) wn = en を解く。
                    wn = torch.linalg.solve(system, eye[:, n])
                except RuntimeError as exc:
                    if "singular" not in str(exc).lower():
                        raise
                    # それでも特異な場合は、擬似逆行列で最小二乗解を得る。
                    wn = torch.linalg.pinv(system) @ eye[:, n]
                # wn^H Vn wn = 1となるよう正規化し、音源の尺度を安定させる。
                norm = torch.sqrt(torch.real(wn.conj() @ VnReg @ wn).clamp_min(threshold))
                wn = wn / norm
                # 分離行列のn行と、その行を適用した推定音源を直ちに更新する。
                W[f, n] = wn.conj()
                Y[f, :, n] = X[f] @ wn.conj()
        if isDraw:
            cost[iIter + 1] = local_calcFdicaCost(Y, W, srcModel)
        if verbose:
            iterationTime = time.perf_counter() - iterationStart
            print(f"[FDICA] 学習 {iIter + 1}/{nIter} 完了 ({iterationTime:.2f}秒)", flush=True)
    return Y, W.permute(1, 2, 0).contiguous(), cost


def local_projectionBack(
    Y: torch.Tensor, S: torch.Tensor, W: torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor]:
    """projection backにより、分離音の振幅を基準マイクの尺度へ戻す。

    ICAでは音源の大きさが一意に決まらない。この処理は、推定音源を混ぜると
    指定した基準マイクの観測信号へ近づくような係数を求め、聞きやすい振幅と
    位相へ補正する。
    """
    if S.ndim == 2:
        S = S[:, :, None]
    nFreq, _, nSrc = Y.shape
    nRef = S.shape[2]
    fixY = torch.empty((nFreq, Y.shape[1], nSrc, nRef), dtype=Y.dtype, device=Y.device)
    fixW = torch.empty((nSrc, W.shape[1], nFreq, nRef), dtype=W.dtype, device=W.device)
    for f in range(nFreq):
        Yf = Y[f].mT  # (source, frame)
        # A = S Y^H (Y Y^H)^-1。pinvはランク落ち時にも計算を継続できる。
        A = S[f].mT @ Yf.mH @ torch.linalg.pinv(Yf @ Yf.mH)  # reference,source
        fixY[f] = (A[:, :, None] * Yf[None, :, :]).permute(2, 1, 0)
        fixW[:, :, f, :] = (A[:, :, None] * W[:, :, f][None, :, :]).permute(1, 2, 0)
    return fixY.squeeze(-1) if nRef == 1 else fixY, fixW.squeeze(-1) if nRef == 1 else fixW


def _time_domain_filter(obsSigInput: torch.Tensor, demixMat: torch.Tensor, fftSize: int) -> torch.Tensor:
    """周波数領域の分離行列をFIRフィルタへ変換し、時間領域で畳み込む。

    ``isFilt=True`` のときに使う。通常の逆STFTとは異なり、各マイク信号へ
    分離フィルタを線形畳み込みし、その和を各推定音源とする。
    """
    # rfft側にない負周波数成分を共役対称性から作り、全周波数へ拡張する。
    full = torch.cat((demixMat, torch.conj(torch.flip(demixMat[:, :, 1:-1], dims=(2,)))), dim=2)
    # 周波数応答を逆FFTすると、時間領域のFIRフィルタ係数になる。
    filt = torch.fft.ifft(full, n=fftSize, dim=2).real
    filt = torch.roll(filt, shifts=fftSize // 2 + 1, dims=2)
    outputs = []
    for n in range(filt.shape[0]):
        # 音源nについて、各マイクの畳み込み結果を足し合わせる。
        total = None
        for ch in range(filt.shape[1]):
            value = torch_functional.conv1d(
                obsSigInput[:, ch][None, None, :],
                torch.flip(filt[n, ch], dims=(0,))[None, None, :],
                padding=fftSize - 1,
            )[0, 0]
            total = value if total is None else total + value
        outputs.append(total)
    return torch.stack(outputs, dim=1)[fftSize // 2 + 1 :]


def bssAuxFdica(
    obsSig: torch.Tensor,
    nSrc: int,
    *,
    fftSize: int = 1024,
    shiftSize: int = 512,
    nIter: int = 50,
    isWhiten: bool = True,
    srcModel: str = "LAP",
    refMic: int | Sequence[int] = 1,
    permSolver: str = "COR",
    isDraw: bool = False,
    sampFreq: float = 16000,
    isPowRatio: bool = True,
    typeCor: str = "Gl+Lo",
    deltaFreq: int = 3,
    ratioFreq: int = 2,
    micPos: torch.Tensor | Sequence[float] | None = None,
    srcSig: torch.Tensor | None = None,
    isFilt: bool = False,
    seed: int = 0,
    verbose: bool = False,
) -> tuple[torch.Tensor, torch.Tensor]:
    """複数マイクの混合信号をMATLAB互換のAuxFDICAで音源分離する。

    引数:
        obsSig: ``(sample, channel)`` 形状の実数観測信号。
        nSrc: 分離したい音源数。マイク数以下の正の整数を指定する。
        fftSize: STFTのFFTサイズ。大きいほど周波数分解能が高くなる。
        shiftSize: STFTのフレームシフト。小さいほど時間方向の重なりが増える。
        nIter: AuxFDICAの反復回数。
        isWhiten: Trueなら周波数ごとにPCA白色化を行う。
        srcModel: ``"LAP"`` または ``"TVG"``。更新に使う音源モデル。
        refMic: projection backの基準マイク番号。MATLABに合わせて1始まり。
        permSolver: ``"none"``, ``"COR"``, ``"DOA"``, ``"IPS"`` のいずれか。
        isDraw: Trueなら各段階の波形・スペクトログラムとコスト履歴を描く。
        sampFreq: サンプリング周波数。描画とDOA計算に使用する。
        isPowRatio: CORソルバーでパワー比を使うかどうか。
        typeCor: CORソルバーの相関種類。``"Gl"``, ``"Lo"``, ``"Gl+Lo"``。
        deltaFreq: CORソルバーが参照する前後の周波数ビン数。
        ratioFreq: CORソルバーが参照する倍音の最大倍率。
        micPos: DOAソルバー用の2本のマイク位置。
        srcSig: IPSソルバー用の正解信号 ``(sample, channel, source)``。
        isFilt: Trueなら逆STFTの代わりに時間領域分離フィルタを使う。
        seed: 再現性確保のための乱数シード。
        verbose: Trueなら処理段階と反復進捗をターミナルへ即時表示する。

    戻り値:
        ``estSig`` と ``cost`` の組。``estSig`` は
        ``(sample, source)`` 形状の分離時間波形。``cost`` は描画時には
        初期値を含む ``nIter + 1`` 個、それ以外ではMATLAB互換の領域を返す。
    """
    # 早い段階で入力を検査し、行列演算の奥で分かりにくいエラーになるのを防ぐ。
    if not isinstance(obsSig, torch.Tensor) or obsSig.ndim != 2 or obsSig.shape[0] == 0:
        raise ValueError("obsSig must be a non-empty tensor shaped (sample, channel)")
    if obsSig.is_complex() or not obsSig.dtype.is_floating_point:
        raise TypeError("obsSig must be a real floating-point tensor")
    if not isinstance(nSrc, int) or nSrc <= 0 or nSrc > obsSig.shape[1]:
        raise ValueError("nSrc must be a positive integer no greater than channel count")
    if not isinstance(nIter, int) or nIter <= 0:
        raise ValueError("nIter must be a positive integer")
    if srcModel not in {"LAP", "TVG"}:
        raise ValueError('srcModel must be "LAP" or "TVG"')
    if permSolver not in {"none", "COR", "DOA", "IPS"}:
        raise ValueError('permSolver must be "none", "COR", "DOA", or "IPS"')
    # 単一の整数も共通処理できるよう、一度リストへそろえる。
    refs = [refMic] if isinstance(refMic, int) else list(refMic)
    if not refs or any(not isinstance(r, int) or r < 1 or r > obsSig.shape[1] for r in refs):
        raise ValueError("refMic contains an invalid one-based microphone index")
    if len(refs) != 1:
        raise ValueError("Python time-domain output currently requires one reference microphone")

    # CPU/GPUで乱数を使う処理が加わっても結果を再現できるようにする。
    torch.manual_seed(seed)
    startTime = time.perf_counter()
    sigLen = obsSig.shape[0]
    if verbose:
        print(
            f"[FDICA] 実行開始: samples={sigLen}, channels={obsSig.shape[1]}, "
            f"sources={nSrc}, device={obsSig.device}",
            flush=True,
        )
        print(f"[FDICA] STFT開始: fftSize={fftSize}, shiftSize={shiftSize}", flush=True)
    # 1. 時間波形を、FDICAが処理する複素スペクトログラムへ変換する。
    obsSpec = dgt_stft(obsSig, fftSize, shiftSize)
    if verbose:
        print(f"[FDICA] STFT完了: spectrum_shape={tuple(obsSpec.shape)}", flush=True)
        print("[FDICA] 白色化開始" if isWhiten else "[FDICA] 白色化をスキップ", flush=True)
    # 2. 必要に応じて白色化し、音源数と同じ次元へそろえる。
    obsSpecInput = local_whitening(obsSpec, nSrc) if isWhiten else obsSpec[:, :, :nSrc]
    if verbose and isWhiten:
        print("[FDICA] 白色化完了", flush=True)
    if verbose:
        print(f"[FDICA] 学習開始: model={srcModel}, iterations={nIter}", flush=True)
    # 3. 補助関数法で周波数ごとの分離行列を反復更新する。
    estSpecFdica, demixMat, cost = local_auxFdica(obsSpecInput, nIter, srcModel, isDraw, verbose)
    if verbose:
        print("[FDICA] 学習完了", flush=True)
        print("[FDICA] Projection back開始", flush=True)
    # 4. ICAで不定となる音量と位相を、指定した基準マイクへ合わせる。
    fixed, demixFixed = local_projectionBack(estSpecFdica, obsSpec[:, :, refs[0] - 1], demixMat)
    if verbose:
        print("[FDICA] Projection back完了", flush=True)
        print(f"[FDICA] 置換問題の解決開始: solver={permSolver}", flush=True)

    # 5. 周波数ごとに入れ替わった音源番号を、選択した方法で統一する。
    if permSolver == "none":
        estSpec = fixed
        estPerm = torch.arange(nSrc, device=obsSig.device).repeat(fixed.shape[0], 1)
    elif permSolver == "COR":
        estSpec, estPerm = permSolverCor(fixed, isPowRatio, typeCor, deltaFreq, ratioFreq)
    elif permSolver == "DOA":
        if micPos is None:
            raise ValueError("micPos is required for the DOA permutation solver")
        estSpec, estPerm = permSolverDoa(demixFixed, fixed, torch.as_tensor(micPos), sampFreq, seed=seed)
    else:
        if srcSig is None or srcSig.ndim != 3:
            raise ValueError("srcSig=(sample, channel, source) is required for IPS")
        srcSpec = dgt_stft(srcSig[:, refs[0] - 1, :], fftSize, shiftSize)
        estSpec, estPerm = permSolverIps(fixed, srcSpec)
    if verbose:
        print("[FDICA] 置換問題の解決完了", flush=True)
    # スペクトログラムだけでなく分離行列にも同じ置換を適用する。
    demixFixed = torch.stack([demixFixed[estPerm[f], :, f] for f in range(fixed.shape[0])], dim=2)

    if verbose:
        method = "時間領域分離フィルタ" if isFilt else "逆STFT"
        print(f"[FDICA] 時間領域信号の生成開始: method={method}", flush=True)
    # 6. 周波数領域の分離結果を、再生可能な時間波形へ戻す。
    if isFilt:
        obsSigInput = dgt_istft(obsSpecInput, fftSize, shiftSize)
        estSig = _time_domain_filter(obsSigInput, demixFixed, fftSize)[:sigLen]
    else:
        estSig = dgt_istft(estSpec, fftSize, shiftSize, length=sigLen)
    if verbose:
        print(f"[FDICA] 時間領域信号の生成完了: output_shape={tuple(estSig.shape)}", flush=True)
    # 7. MATLAB版と同じ処理段階を可視化する。showは全ウィンドウを閉じるまで待つ。
    if isDraw:
        if verbose:
            print("[FDICA] 描画開始: ウィンドウを閉じるまで処理を待機します", flush=True)
        local_plotSpectrogram(obsSig, sampFreq, fftSize, shiftSize, title="Observed signal")
        local_plotSpectrogram(dgt_istft(obsSpecInput, fftSize, shiftSize, length=sigLen), sampFreq, fftSize, shiftSize, title="FDICA input signal")
        local_plotSpectrogram(dgt_istft(estSpecFdica, fftSize, shiftSize, length=sigLen), sampFreq, fftSize, shiftSize, title="Estimated signal before projection back")
        if permSolver != "none":
            local_plotSpectrogram(dgt_istft(fixed, fftSize, shiftSize, length=sigLen), sampFreq, fftSize, shiftSize, title="Estimated signal before permutation solver")
        local_plotSpectrogram(estSig, sampFreq, fftSize, shiftSize, title="Estimated signal")
        local_plotCost(cost, nIter)
        local_showPlots()
        if verbose:
            print("[FDICA] 描画完了", flush=True)
    if verbose:
        elapsed = time.perf_counter() - startTime
        print(f"[FDICA] 実行完了 ({elapsed:.2f}秒)", flush=True)
    return estSig, cost
