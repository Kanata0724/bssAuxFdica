# 対応注: 学習用コピー：元ファイル python_fdica/bssAuxFdica.py（learning / b934fef 時点）。
# 対応注: import、型注釈、device/dtype、時間測定、verbose表示はPython固有。描画は機能単位で対応。
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


# [bssAuxFdica.m:190] function [Y, dP] = local_whitening(X, N)
# 対応注: PythonはYのみ返す。
def local_whitening(X: torch.Tensor, nSrc: int) -> torch.Tensor:
    """周波数ごとにPCA白色化し、成分間の相関とスケール差を取り除く。

    入力 ``X`` は ``(frequency, frame, channel)``、出力 ``Y`` は
    ``(frequency, frame, source)``。白色化後の共分散行列は単位行列に
    近くなり、FDICAが分離行列を求めやすくなる。
    """
    # [bssAuxFdica.m:202] [I, J, ~] = size(X, [1, 2, 3]); % nFreq x nTime x nCh
    # 対応注: I=nFreq、J=nFrame、N=nSrc。
    nFreq, nFrame, _ = X.shape
    # [bssAuxFdica.m:203] Y = zeros(I, J, N);
    # 対応注: Pythonは未初期化領域を確保し、各周波数で全要素を書き込む。
    Y = torch.empty((nFreq, nFrame, nSrc), dtype=X.dtype, device=X.device)
    # 対応注: Python独自：固有値の下限値。MATLAB版にはこの下限制御がない。
    eps = torch.finfo(X.real.dtype).eps
    # [bssAuxFdica.m:207] for i = 1:I
    # 対応注: MATLABのi=1..Iに対しPythonは0..I-1。
    for iFreq in range(nFreq):
        # 1周波数分を、線形代数で一般的な (channel, frame) へ転置する。
        # [bssAuxFdica.m:208] Xi = Xp(:, :, i); % M x J
        # 対応注: MATLABでは先にXp=permute(X,[3,2,1])。Pythonは周波数を選んでから非共役転置。
        Xi = X[iFreq].mT
        # Xi Xi^H / フレーム数は、チャンネル間の複素共分散行列である。
        # [bssAuxFdica.m:209] V = Xi*(Xi')/J; % covariance matrix of data matrix X (K x K)
        # 対応注: 平均を引かない二次モーメント。Xi.mHは複素共役転置。
        covariance = Xi @ Xi.mH / nFrame
        # Hermite行列専用のeighを使い、固有値と直交固有ベクトルを求める。
        # [bssAuxFdica.m:210] [P, D] = eig(V); % eigenvalue decomposition (V = P*D*inv(P), P includes eigenvectors and D is a diagonal matrix with eigenvalues)
        # 対応注: PythonはHermite行列用eigh。固有値は対角行列でなくベクトル。
        eigenvalues, eigenvectors = torch.linalg.eigh(covariance)
        # 固有値が大きい順にnSrc成分を残す。これは主成分分析に相当する。
        # [bssAuxFdica.m:211] [~, idx] = sort(diag(D), "descend"); % sort eigenvalues in descending order
        # [bssAuxFdica.m:212] D = D(idx, idx); % sorted D
        # [bssAuxFdica.m:213] P = P(:, idx); % sorted P
        # 対応注: Pythonはこの時点で上位nSrc個だけ選ぶ。
        idx = torch.argsort(eigenvalues, descending=True)[:nSrc]
        # [bssAuxFdica.m:214] dP = P(:, 1:N); % top-d eigenvectors
        dP = eigenvectors[:, idx]
        # 主成分へ射影し、各成分を固有値の平方根で割って分散を1へそろえる。
        # [bssAuxFdica.m:215] Yi = sqrt(D)\(dP')*Xi; % whitened vector (N x J)
        # [bssAuxFdica.m:216] Y(i, :, :) = Yi.'; % J x N
        # 対応注: 処理目的の対応。Pythonは選択固有値のみを使いepsで下限を設ける。MATLAB原文はDをN次元に切り詰めず、N<Mの場合に次元不整合となる。
        Y[iFreq] = ((dP.mH @ Xi) / torch.sqrt(eigenvalues[idx].clamp_min(eps))[:, None]).mT
    return Y


# [bssAuxFdica.m:290] function costVal = local_calcFdicaCost(Yp, W, srcModel, I, J)
def local_calcFdicaCost(Y: torch.Tensor, W: torch.Tensor, srcModel: str) -> torch.Tensor:
    """FDICAの目的関数（コスト関数）を計算する。

    ``Y`` は ``(frequency, frame, source)``、``W`` は周波数ごとの分離行列。
    音源モデルに応じた罰則と、分離行列が退化しないためのlog行列式項から
    構成される。値そのものより、反復に伴う変化を見るために利用する。
    """
    # [bssAuxFdica.m:290] function costVal = local_calcFdicaCost(Yp, W, srcModel, I, J)
    # 対応注: MATLAB引数JをPythonではY.shape[1]から取得。
    nFrame = Y.shape[1]
    # log(0)を避けるため、行列式の絶対値を機械イプシロン以上に制限する。
    # [bssAuxFdica.m:293] detW(i,1) = det(W(:,:,i));
    # 対応注: 全周波数のdetとlog(abs(det))の総和をまとめて計算。Pythonだけdet絶対値にeps下限を設定。
    logdet = torch.log(torch.abs(torch.linalg.det(W)).clamp_min(torch.finfo(Y.real.dtype).eps)).sum()
    # [bssAuxFdica.m:295] if srcModel == "LAP"
    # 対応注: LAP分岐（以下のreturnがMATLABのcostVal代入に対応）。
    if srcModel == "LAP":
        # [bssAuxFdica.m:296] costVal = sum(abs(Yp), "all") - 2*J*sum(log(abs(detW)));
        return torch.abs(Y).sum() - 2 * nFrame * logdet
    # [bssAuxFdica.m:298] costVal = sum(log(max(abs(Yp).^2, eps)), "all") - 2*J*sum(log(abs(detW)));
    return torch.log((torch.abs(Y) ** 2).clamp_min(torch.finfo(Y.real.dtype).eps)).sum() - 2 * nFrame * logdet


# [DGTtool.m:567] function plot(obj,x,fs,options)
# 対応注: 機能単位の対応。MATLABの描画をmatplotlibで再構成しており、各行の一対一対応ではない。
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
    # [bssAuxFdica.m:359] figure; plot(0:nIter, cost);
    # 対応注: 次行のax.plotまでが対応。
    fig, ax = plt.subplots()
    ax.plot(range(nIter + 1), costCpu[: nIter + 1].numpy())
    # [bssAuxFdica.m:361] xlabel("Number of iterations"); ylabel("Value of cost function");
    # 対応注: Pythonではxlabelとylabelを別々に呼ぶ。
    ax.set_xlabel("Number of iterations")
    ax.set_ylabel("Value of cost function")
    # [bssAuxFdica.m:362] grid on;
    ax.grid(True)
    ax.tick_params(labelsize=12)


# 対応注: Python独自：matplotlibのウィンドウを表示して閉じるまで待機。
def local_showPlots() -> None:
    """作成した全Figureを表示し、利用者が閉じるまで処理を待つ。"""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("isDraw=True requires matplotlib. Install it with `pip install matplotlib`.") from exc

    plt.show(block=True)


# [bssAuxFdica.m:221] function [Y, W, cost] = local_auxFdica(X, nIter, srcModel, isDraw)
# 対応注: 内部WはPythonが周波数×音源×channel、MATLABが音源×channel×周波数。
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
    # [bssAuxFdica.m:236] [I, J, M] = size(X, [1,2,3]); % nFreq x nTime x nCh
    # [bssAuxFdica.m:237] N = M; % number of sources
    # 対応注: I=nFreq、J=nFrame、M=N=nCh。
    nFreq, nFrame, nCh = X.shape
    # 最初は「入力をそのまま出力する」単位行列を全周波数へ用意する。
    # [bssAuxFdica.m:238] E = repmat(eye(M), [1, 1, I]);
    # [bssAuxFdica.m:239] W = E; % initial demixing matrix (N x M x I)
    # 対応注: Pythonは周波数軸を先頭に置く。
    W = torch.eye(nCh, dtype=X.dtype, device=X.device).repeat(nFreq, 1, 1)
    # [bssAuxFdica.m:240] Y = X; % initial estimated spectrogram
    # 対応注: PythonのYはI×J×Nのまま。MATLABは後でYp=N×J×Iに並べ替える。
    Y = X.clone()
    # MATLAB版と同様、描画するときだけ初期値を含むコストを計算する。
    # [bssAuxFdica.m:244] cost = zeros(nIter, 1);
    # 対応注: Pythonは描画時のnIter+1要素を先に確保。MATLABは代入時に配列を拡張する。
    cost = torch.zeros(nIter + 1 if isDraw else nIter, dtype=X.real.dtype, device=X.device)
    if isDraw:
        # [bssAuxFdica.m:246] cost(1,1) = local_calcFdicaCost(Yp, W, srcModel, I, J);
        cost[0] = local_calcFdicaCost(Y, W, srcModel)
    # 0除算や極端に大きな更新を避けるための、dtypeに応じた下限値。
    # [bssAuxFdica.m:254] Rp = max(abs(Yp), 10000*eps);
    # 対応注: この行の下限10000*epsをPythonでは先に変数化。
    threshold = 10000 * torch.finfo(X.real.dtype).eps
    # [bssAuxFdica.m:238] E = repmat(eye(M), [1, 1, I]);
    # 対応注: Pythonは全周波数で共通の単位行列を使い回す。
    eye = torch.eye(nCh, dtype=X.dtype, device=X.device)
    # [bssAuxFdica.m:251] for iIter = 1:nIter
    for iIter in range(nIter):
        # 対応注: Python独自：処理時間の測定。verboseによる進捗表示もMATLABのfprintfと表示仕様が異なる。
        iterationStart = time.perf_counter()
        if verbose:
            print(f"[FDICA] 学習 {iIter + 1}/{nIter} 開始", flush=True)
        # LAPは振幅、TVGはパワーを補助変数の重みに使う。
        # [bssAuxFdica.m:254] Rp = max(abs(Yp), 10000*eps);
        # [bssAuxFdica.m:255] elseif srcModel == "TVG"
        # [bssAuxFdica.m:256] Rp = max(abs(Yp).^2, 10000*eps);
        # 対応注: LAPは絶対値、TVGは絶対値の二乗。下限処理は次行のclamp_minで行う。
        radius = torch.abs(Y) if srcModel == "LAP" else torch.abs(Y) ** 2
        # [bssAuxFdica.m:259] invRp = 1./Rp; % N x J x I
        # 対応注: Rpの下限制御を含む。重みは外側反復ごとに計算し、音源更新中は再計算しない。
        invRadius = radius.clamp_min(threshold).reciprocal()
        # [bssAuxFdica.m:260] for n = 1:N
        for n in range(nCh):
            # [bssAuxFdica.m:262] Vk = pagemtimes(D.*Xp, Xph)/J; % M x M x I, pagewise matrix multiplication ((D(:,:,i).*Xp(:,:,i))*Xp(:,:,i)'/J)
            # 対応注: MATLABのページ演算を、Pythonは周波数ループ内で実行する。
            for f in range(nFreq):
                # 音源n・周波数fに対応する重み付き共分散行列Vnを作る。
                # [bssAuxFdica.m:241] Xp = permute(X, [3, 2, 1]); % M x J x I
                # 対応注: 周波数fだけ取り出してchannel×frameにする。Xf.mHがMATLABのXphの各ページに対応。
                Xf = X[f].mT  # (channel, frame)
                # [bssAuxFdica.m:261] D = repmat(invRp(n, :, :), [M, 1, 1]); % M x J x I
                # [bssAuxFdica.m:262] Vk = pagemtimes(D.*Xp, Xph)/J; % M x M x I, pagewise matrix multiplication ((D(:,:,i).*Xp(:,:,i))*Xp(:,:,i)'/J)
                # 対応注: Pythonは重みの複製をbroadcastingで代替。
                Vn = (Xf * invRadius[f, :, n][None, :]) @ Xf.mH / nFrame
                # 無音に近い帯域や強く相関したマイク信号ではVnが特異行列に
                # なり得る。微小値を対角へ加え、数値的に解ける状態へする。
                # 対応注: Python独自：Vnに微小な対角成分ridgeを追加する。MATLAB本体に対応する追加処理はない。
                meanPower = torch.real(torch.trace(Vn)) / nCh
                ridge = 100 * torch.finfo(X.real.dtype).eps * meanPower.clamp_min(1.0)
                VnReg = Vn + ridge * eye
                # [bssAuxFdica.m:263] wn = pagemldivide(pagemtimes(W, Vk), E(:, n, :)); % M x 1 x I, pagewise operation ((W(:,:,i)*Vk(:,:,i)) \ E(:, n, :))
                # 対応注: 以下のsolveと合わせて左除算に対応。ただしPythonはVkではなく正則化したVnRegを使用。
                system = W[f] @ VnReg
                try:
                    # iterative projectionの更新式 (W Vn) wn = en を解く。
                    wn = torch.linalg.solve(system, eye[:, n])
                # 対応注: Python独自：特異行列の場合だけpinvによる解へ切り替える。MATLAB本体に明示的な同じ分岐はない。
                except RuntimeError as exc:
                    if "singular" not in str(exc).lower():
                        raise
                    # それでも特異な場合は、擬似逆行列で最小二乗解を得る。
                    wn = torch.linalg.pinv(system) @ eye[:, n]
                # wn^H Vn wn = 1となるよう正規化し、音源の尺度を安定させる。
                # [bssAuxFdica.m:264] wn = wn ./ sqrt( pagemtimes(pagemtimes(wn, "ctranspose", Vk, "none"), wn) ); % M x 1 x I, pagewise operation (wn(:,:,i)/sqrt(wn(:,:,i)'*Vk(:,:,i)*wn(:,:,i)))
                # 対応注: 次行の除算までで対応。PythonはVnReg、実部抽出、下限制御を使うため数値処理は同一ではない。
                norm = torch.sqrt(torch.real(wn.conj() @ VnReg @ wn).clamp_min(threshold))
                wn = wn / norm
                # 分離行列のn行と、その行を適用した推定音源を直ちに更新する。
                # [bssAuxFdica.m:267] W(n, :, :) = wnh;
                # 対応注: wnh=wnの共役転置。Pythonは一次元ベクトルを行へ代入。MATLABではYp更新の後にWを代入。
                W[f, n] = wn.conj()
                # [bssAuxFdica.m:266] Yp(n, :, :) = pagemtimes(wnh, Xp); % 1 x J x I, pagewise matrix multiplication (wnh(:,:,i)*Xp(:,:,i))
                # 対応注: 両者ともwn^H xを計算。保持する配列の軸順が違う。
                Y[f, :, n] = X[f] @ wn.conj()
        if isDraw:
            # [bssAuxFdica.m:282] cost(iIter+1, 1) = local_calcFdicaCost(Yp, W, srcModel, I, J);
            # 対応注: 反復番号はPythonが0始まりなので、同じ反復の保存先は添字表現が異なる。
            cost[iIter + 1] = local_calcFdicaCost(Y, W, srcModel)
        if verbose:
            iterationTime = time.perf_counter() - iterationStart
            print(f"[FDICA] 学習 {iIter + 1}/{nIter} 完了 ({iterationTime:.2f}秒)", flush=True)
    # [bssAuxFdica.m:285] Y = permute(Yp, [3, 2, 1]);
    # 対応注: 出力形式をそろえる処理。PythonのYは既にI×J×Nで、Wの軸を戻す。
    return Y, W.permute(1, 2, 0).contiguous(), cost


# [bssAuxFdica.m:303] function [fixY, fixW] = local_projectionBack(Y, S, W)
def local_projectionBack(
    Y: torch.Tensor, S: torch.Tensor, W: torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor]:
    """projection backにより、分離音の振幅を基準マイクの尺度へ戻す。

    ICAでは音源の大きさが一意に決まらない。この処理は、推定音源を混ぜると
    指定した基準マイクの観測信号へ近づくような係数を求め、聞きやすい振幅と
    位相へ補正する。
    """
    # 対応注: Pythonの二次元参照信号を三次元にする。MATLABでは末尾の長さ1の次元を暗黙に扱える。
    if S.ndim == 2:
        S = S[:, :, None]
    nFreq, _, nSrc = Y.shape
    nRef = S.shape[2]
    fixY = torch.empty((nFreq, Y.shape[1], nSrc, nRef), dtype=Y.dtype, device=Y.device)
    fixW = torch.empty((nSrc, W.shape[1], nFreq, nRef), dtype=W.dtype, device=W.device)
    for f in range(nFreq):
        # [bssAuxFdica.m:321] Yp = permute(Y, [3, 2, 1]); % N x J x I
        # 対応注: Pythonは周波数ごとにsource×frameへ変換。
        Yf = Y[f].mT  # (source, frame)
        # A = S Y^H (Y Y^H)^-1。pinvはランク落ち時にも計算を継続できる。
        # [bssAuxFdica.m:324] Yph = pagectranspose(Yp); % J x N x I, pagewise Hermitian transpose (Yp')
        # [bssAuxFdica.m:325] YpYph = pagemtimes(Yp, Yph); % N x N x I, pagewise matrix multiplication (Yp*Yp')
        # [bssAuxFdica.m:326] YphOnYpYph = pagemrdivide(Yph, YpYph); % J x N x I, pagewise matrix right-division (Yp'/(Yp*Yp'))
        # [bssAuxFdica.m:327] A = pagemtimes(Sp, YphOnYpYph); % 1 x N x I or M x N x I, pagewise matrix multiplication (Sp * Yp'/(Yp*Yp'))
        # 対応注: 処理ブロックの対応。MATLABの右除算に対しPythonはGram行列のpinvを使用。ランク落ち時まで同値とは限らない。
        A = S[f].mT @ Yf.mH @ torch.linalg.pinv(Yf @ Yf.mH)  # reference,source
        # [bssAuxFdica.m:330] fixY = Ap .* Ypp; % M x N x J x I, using implicit expansion
        # [bssAuxFdica.m:331] fixY = permute(fixY, [4, 3, 2, 1]); % I x J x N x M
        # 対応注: 係数Aを音源・参照マイク別に乗算し、軸を戻す。
        fixY[f] = (A[:, :, None] * Yf[None, :, :]).permute(2, 1, 0)
        # [bssAuxFdica.m:332] fixW = Ap .* Wp; % M x N x N x I, using implicit expansion
        # [bssAuxFdica.m:333] fixW = permute(fixW, [2, 3, 4, 1]); % N x N x I x M
        # 対応注: Wの各音源行にも同じ係数を乗算。
        fixW[:, :, f, :] = (A[:, :, None] * W[:, :, f][None, :, :]).permute(1, 2, 0)
    # 対応注: Pythonでは参照マイクが1本のとき末尾軸を除去。MATLABは末尾の長さ1の次元を暗黙に扱う。
    return fixY.squeeze(-1) if nRef == 1 else fixY, fixW.squeeze(-1) if nRef == 1 else fixW


# [bssAuxFdica.m:155] if isFilt
# 対応注: MATLAB本体のisFilt分岐をPythonでは別関数に抽出。
def _time_domain_filter(obsSigInput: torch.Tensor, demixMat: torch.Tensor, fftSize: int) -> torch.Tensor:
    """周波数領域の分離行列をFIRフィルタへ変換し、時間領域で畳み込む。

    ``isFilt=True`` のときに使う。通常の逆STFTとは異なり、各マイク信号へ
    分離フィルタを線形畳み込みし、その和を各推定音源とする。
    """
    # rfft側にない負周波数成分を共役対称性から作り、全周波数へ拡張する。
    # [bssAuxFdica.m:158] W = cat(3, demixMatFix, flip(conj(demixMatFix(:, :, 2:end-1)), 3)); % produce beyond Nyquist components
    full = torch.cat((demixMat, torch.conj(torch.flip(demixMat[:, :, 1:-1], dims=(2,)))), dim=2)
    # 周波数応答を逆FFTすると、時間領域のFIRフィルタ係数になる。
    # [bssAuxFdica.m:159] demixFilt = real(ifft(W, fftSize, 3)); % fftSize x nSrc x nMic
    filt = torch.fft.ifft(full, n=fftSize, dim=2).real
    # [bssAuxFdica.m:160] demixFilt = circshift(demixFilt, fftSize/2+1, 3); % move peak to center by circular shifting
    filt = torch.roll(filt, shifts=fftSize // 2 + 1, dims=2)
    outputs = []
    # [bssAuxFdica.m:161] for iSrc = 1:nSrc
    for n in range(filt.shape[0]):
        # 音源nについて、各マイクの畳み込み結果を足し合わせる。
        total = None
        # [bssAuxFdica.m:162] for iCh = 1:nCh
        # 対応注: Pythonはフィルタのchannel数、MATLABは元観測のnChを使用。次元削減時は要注意。
        for ch in range(filt.shape[1]):
            # [bssAuxFdica.m:163] f = squeeze(demixFilt(iSrc, iCh, :));
            # [bssAuxFdica.m:164] tmp(:, iCh) = conv(obsSigInput(:, iCh), f); % linear convolution
            # 対応注: conv1dは相互相関なので係数をflipし、paddingで完全な線形畳み込みを得る。
            value = torch_functional.conv1d(
                obsSigInput[:, ch][None, None, :],
                torch.flip(filt[n, ch], dims=(0,))[None, None, :],
                padding=fftSize - 1,
            )[0, 0]
            # [bssAuxFdica.m:166] estSig(:, iSrc) = sum(tmp, 2);
            # 対応注: Pythonはchannelループ内で逐次加算。
            total = value if total is None else total + value
        outputs.append(total)
    # [bssAuxFdica.m:168] estSig(1:fftSize/2+1,:) = []; % cut initial components caused by group delay (circular shifting)
    # 対応注: 音源を列に並べ、先頭の遅延成分を除く。
    return torch.stack(outputs, dim=1)[fftSize // 2 + 1 :]


# [bssAuxFdica.m:1] function [estSig, cost] = bssAuxFdica(obsSig, nSrc, args)
# 対応注: 名前付き引数はMATLABのargumentsブロックに対応。seedとverboseはPython独自の関数引数。
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
    # [bssAuxFdica.m:85] obsSig (:,:) double
    # [bssAuxFdica.m:86] nSrc (1,1) double {mustBeInteger, mustBePositive}
    # [bssAuxFdica.m:87] args.fftSize (1,1) double {mustBeInteger, mustBePositive} = 1024
    # 対応注: 以下はMATLABのargumentsとCheck argument errorsに対応する入力検査。dtypeなど検査条件は完全には一致しない。
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
    # 対応注: refMicは両版とも外部指定は1始まり。Pythonでは参照時に1を引く。Python公開関数は参照マイク1本のみ許可。
    refs = [refMic] if isinstance(refMic, int) else list(refMic)
    if not refs or any(not isinstance(r, int) or r < 1 or r > obsSig.shape[1] for r in refs):
        raise ValueError("refMic contains an invalid one-based microphone index")
    if len(refs) != 1:
        raise ValueError("Python time-domain output currently requires one reference microphone")

    # CPU/GPUで乱数を使う処理が加わっても結果を再現できるようにする。
    # [main.m:25] rng(seed);
    # 対応注: MATLABではmain側、Pythonでは関数内部でseedを設定。乱数列そのものの一致は保証しない。
    torch.manual_seed(seed)
    startTime = time.perf_counter()
    # [bssAuxFdica.m:115] [sigLen, nCh] = size(obsSig, [1, 2]);
    sigLen = obsSig.shape[0]
    if verbose:
        print(
            f"[FDICA] 実行開始: samples={sigLen}, channels={obsSig.shape[1]}, "
            f"sources={nSrc}, device={obsSig.device}",
            flush=True,
        )
        print(f"[FDICA] STFT開始: fftSize={fftSize}, shiftSize={shiftSize}", flush=True)
    # 1. 時間波形を、FDICAが処理する複素スペクトログラムへ変換する。
    # [bssAuxFdica.m:121] F = DGTtool("windowName", "b", "windowLength", fftSize, "windowShift", shiftSize); % create DGTtool instance
    # [bssAuxFdica.m:122] obsSpec = F.DGT(obsSig); % STFT
    # 対応注: PythonではDGTtoolオブジェクトを作らず関数を呼ぶ。
    obsSpec = dgt_stft(obsSig, fftSize, shiftSize)
    if verbose:
        print(f"[FDICA] STFT完了: spectrum_shape={tuple(obsSpec.shape)}", flush=True)
        print("[FDICA] 白色化開始" if isWhiten else "[FDICA] 白色化をスキップ", flush=True)
    # 2. 必要に応じて白色化し、音源数と同じ次元へそろえる。
    # [bssAuxFdica.m:126] if isWhiten
    # [bssAuxFdica.m:127] obsSpecInput = local_whitening(obsSpec, nSrc);
    # [bssAuxFdica.m:128] else
    # [bssAuxFdica.m:129] obsSpecInput = obsSpec(:, :, 1:nSrc); % discard unnecessary channels
    # [bssAuxFdica.m:130] end
    # 対応注: 1:nSrcはPythonでは:nSrc。
    obsSpecInput = local_whitening(obsSpec, nSrc) if isWhiten else obsSpec[:, :, :nSrc]
    if verbose and isWhiten:
        print("[FDICA] 白色化完了", flush=True)
    if verbose:
        print(f"[FDICA] 学習開始: model={srcModel}, iterations={nIter}", flush=True)
    # 3. 補助関数法で周波数ごとの分離行列を反復更新する。
    # [bssAuxFdica.m:133] [estSpecFdica, demixMat, cost] = local_auxFdica(obsSpecInput, nIter, srcModel, isDraw);
    estSpecFdica, demixMat, cost = local_auxFdica(obsSpecInput, nIter, srcModel, isDraw, verbose)
    if verbose:
        print("[FDICA] 学習完了", flush=True)
        print("[FDICA] Projection back開始", flush=True)
    # 4. ICAで不定となる音量と位相を、指定した基準マイクへ合わせる。
    # [bssAuxFdica.m:136] [estSpecFdicaFix, demixMatFix] = local_projectionBack(estSpecFdica, obsSpec(:,:,refMic), demixMat);
    fixed, demixFixed = local_projectionBack(estSpecFdica, obsSpec[:, :, refs[0] - 1], demixMat)
    if verbose:
        print("[FDICA] Projection back完了", flush=True)
        print(f"[FDICA] 置換問題の解決開始: solver={permSolver}", flush=True)

    # 5. 周波数ごとに入れ替わった音源番号を、選択した方法で統一する。
    # [bssAuxFdica.m:139] if permSolver == "none"
    if permSolver == "none":
        # [bssAuxFdica.m:140] estSpec = estSpecFdicaFix;
        estSpec = fixed
        # [bssAuxFdica.m:141] estPerm = repmat(1:nSrc, [nFreq, 1]);
        # 対応注: 置換番号はMATLABが1始まり、Pythonは0始まり。
        estPerm = torch.arange(nSrc, device=obsSig.device).repeat(fixed.shape[0], 1)
    elif permSolver == "COR":
        # [bssAuxFdica.m:143] [estSpec, estPerm] = permSolverCor(estSpecFdicaFix, args.isPowRatio, args.typeCor, args.deltaFreq, args.ratioFreq);
        estSpec, estPerm = permSolverCor(fixed, isPowRatio, typeCor, deltaFreq, ratioFreq)
    elif permSolver == "DOA":
        if micPos is None:
            raise ValueError("micPos is required for the DOA permutation solver")
        # [bssAuxFdica.m:145] [estSpec, estPerm] = permSolverDoa(demixMatFix, estSpecFdicaFix, args.micPos, args.sampFreq);
        # 対応注: DOA関数内部のクラスタ初期値や無効値処理には実装差がある。
        estSpec, estPerm = permSolverDoa(demixFixed, fixed, torch.as_tensor(micPos), sampFreq, seed=seed)
    else:
        if srcSig is None or srcSig.ndim != 3:
            raise ValueError("srcSig=(sample, channel, source) is required for IPS")
        # [bssAuxFdica.m:147] srcSpect = F.DGT(squeeze(args.srcSig(:, args.refMic, :)));
        # 対応注: Pythonの整数添字はmic軸を取り除くので、この箇所にsqueezeは不要。
        srcSpec = dgt_stft(srcSig[:, refs[0] - 1, :], fftSize, shiftSize)
        # [bssAuxFdica.m:148] [estSpec, estPerm] = permSolverIps(estSpecFdicaFix, srcSpect);
        estSpec, estPerm = permSolverIps(fixed, srcSpec)
    if verbose:
        print("[FDICA] 置換問題の解決完了", flush=True)
    # スペクトログラムだけでなく分離行列にも同じ置換を適用する。
    # [bssAuxFdica.m:150] for iFreq = 1:nFreq
    # [bssAuxFdica.m:151] demixMatFix(:, :, iFreq) = demixMatFix(estPerm(iFreq, :), :, iFreq);
    # [bssAuxFdica.m:152] end
    # 対応注: 各周波数で分離行列の音源行を並べ替える。
    demixFixed = torch.stack([demixFixed[estPerm[f], :, f] for f in range(fixed.shape[0])], dim=2)

    if verbose:
        method = "時間領域分離フィルタ" if isFilt else "逆STFT"
        print(f"[FDICA] 時間領域信号の生成開始: method={method}", flush=True)
    # 6. 周波数領域の分離結果を、再生可能な時間波形へ戻す。
    # [bssAuxFdica.m:155] if isFilt
    if isFilt:
        # [bssAuxFdica.m:157] obsSigInput = F.pinv(obsSpecInput); % observed signal
        obsSigInput = dgt_istft(obsSpecInput, fftSize, shiftSize)
        # [bssAuxFdica.m:173] estSig = estSig(1:sigLen, :);
        # 対応注: 時間領域フィルタ処理は別関数へ移し、その後に元の信号長へ切る。
        estSig = _time_domain_filter(obsSigInput, demixFixed, fftSize)[:sigLen]
    else:
        # [bssAuxFdica.m:171] estSig = F.pinv(estSpec);
        # 対応注: length=sigLenでMATLAB後続の長さ切り詰めも行う。
        estSig = dgt_istft(estSpec, fftSize, shiftSize, length=sigLen)
    if verbose:
        print(f"[FDICA] 時間領域信号の生成完了: output_shape={tuple(estSig.shape)}", flush=True)
    # 7. MATLAB版と同じ処理段階を可視化する。showは全ウィンドウを閉じるまで待つ。
    if isDraw:
        if verbose:
            print("[FDICA] 描画開始: ウィンドウを閉じるまで処理を待機します", flush=True)
        # [bssAuxFdica.m:177] F.plot(obsSig, args.sampFreq); % observed signal
        # 対応注: 続く描画群もMATLAB本体の末尾に対応。描画の細部は独自実装。
        local_plotSpectrogram(obsSig, sampFreq, fftSize, shiftSize, title="Observed signal")
        local_plotSpectrogram(dgt_istft(obsSpecInput, fftSize, shiftSize, length=sigLen), sampFreq, fftSize, shiftSize, title="FDICA input signal")
        local_plotSpectrogram(dgt_istft(estSpecFdica, fftSize, shiftSize, length=sigLen), sampFreq, fftSize, shiftSize, title="Estimated signal before projection back")
        if permSolver != "none":
            local_plotSpectrogram(dgt_istft(fixed, fftSize, shiftSize, length=sigLen), sampFreq, fftSize, shiftSize, title="Estimated signal before permutation solver")
        local_plotSpectrogram(estSig, sampFreq, fftSize, shiftSize, title="Estimated signal")
        # [bssAuxFdica.m:184] local_plotCost(cost, nIter); % cost function behavior
        local_plotCost(cost, nIter)
        local_showPlots()
        if verbose:
            print("[FDICA] 描画完了", flush=True)
    if verbose:
        elapsed = time.perf_counter() - startTime
        print(f"[FDICA] 実行完了 ({elapsed:.2f}秒)", flush=True)
    return estSig, cost
