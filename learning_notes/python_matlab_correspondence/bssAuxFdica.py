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
def local_whitening(X: torch.Tensor, nSrc: int) -> torch.Tensor:
    # [bssAuxFdica.m:202] [I, J, ~] = size(X, [1, 2, 3]); % nFreq x nTime x nCh
    nFreq, nFrame, _ = X.shape
    # [bssAuxFdica.m:203] Y = zeros(I, J, N);
    Y = torch.empty((nFreq, nFrame, nSrc), dtype=X.dtype, device=X.device)
    eps = torch.finfo(X.real.dtype).eps
    # [bssAuxFdica.m:207] for i = 1:I
    for iFreq in range(nFreq):
        # [bssAuxFdica.m:208] Xi = Xp(:, :, i); % M x J
        Xi = X[iFreq].mT
        # [bssAuxFdica.m:209] V = Xi*(Xi')/J; % covariance matrix of data matrix X (K x K)
        covariance = Xi @ Xi.mH / nFrame
        # [bssAuxFdica.m:210] [P, D] = eig(V); % eigenvalue decomposition (V = P*D*inv(P), P includes eigenvectors and D is a diagonal matrix with eigenvalues)
        eigenvalues, eigenvectors = torch.linalg.eigh(covariance)
        # [bssAuxFdica.m:211] [~, idx] = sort(diag(D), "descend"); % sort eigenvalues in descending order
        # [bssAuxFdica.m:212] D = D(idx, idx); % sorted D
        # [bssAuxFdica.m:213] P = P(:, idx); % sorted P
        idx = torch.argsort(eigenvalues, descending=True)[:nSrc]
        # [bssAuxFdica.m:214] dP = P(:, 1:N); % top-d eigenvectors
        dP = eigenvectors[:, idx]
        # [bssAuxFdica.m:215] Yi = sqrt(D)\(dP')*Xi; % whitened vector (N x J)
        # [bssAuxFdica.m:216] Y(i, :, :) = Yi.'; % J x N
        Y[iFreq] = ((dP.mH @ Xi) / torch.sqrt(eigenvalues[idx].clamp_min(eps))[:, None]).mT
    return Y


# [bssAuxFdica.m:290] function costVal = local_calcFdicaCost(Yp, W, srcModel, I, J)
def local_calcFdicaCost(Y: torch.Tensor, W: torch.Tensor, srcModel: str) -> torch.Tensor:
    # [bssAuxFdica.m:290] function costVal = local_calcFdicaCost(Yp, W, srcModel, I, J)
    nFrame = Y.shape[1]
    # [bssAuxFdica.m:293] detW(i,1) = det(W(:,:,i));
    logdet = torch.log(torch.abs(torch.linalg.det(W)).clamp_min(torch.finfo(Y.real.dtype).eps)).sum()
    # [bssAuxFdica.m:295] if srcModel == "LAP"
    if srcModel == "LAP":
        # [bssAuxFdica.m:296] costVal = sum(abs(Yp), "all") - 2*J*sum(log(abs(detW)));
        return torch.abs(Y).sum() - 2 * nFrame * logdet
    # [bssAuxFdica.m:298] costVal = sum(log(max(abs(Yp).^2, eps)), "all") - 2*J*sum(log(abs(detW)));
    return torch.log((torch.abs(Y) ** 2).clamp_min(torch.finfo(Y.real.dtype).eps)).sum() - 2 * nFrame * logdet


# [DGTtool.m:567] function plot(obj,x,fs,options)
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
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("isDraw=True requires matplotlib. Install it with `pip install matplotlib`.") from exc

    if signal.ndim != 2:
        raise ValueError("signal must be shaped (sample, channel/source)")
    x = signal.detach().cpu()
    nSample, nSignal = x.shape
    window = torch.blackman_window(fftSize, periodic=True, dtype=x.dtype)
    normConst = torch.sum(window).item() / 2 if normalize else 1.0
    spec = dgt_stft(x, fftSize, shiftSize).detach().cpu() / normConst
    xForFft = x / normConst
    eps = torch.finfo(spec.real.dtype).eps

    if sampFreq == 1:
        fsPlot = 1.0
        freqUnit = "[periods/sample]"
        timeUnit = "[samples]"
        time = torch.arange(spec.shape[1], dtype=x.real.dtype) * shiftSize
        waveTime = torch.arange(nSample, dtype=x.real.dtype)
        xLimit = (0, max(nSample - 1, 1))
    else:
        fsPlot = sampFreq / 1000
        freqUnit = "[kHz]"
        timeUnit = "[s]"
        time = torch.arange(spec.shape[1], dtype=x.real.dtype) * shiftSize / sampFreq
        waveTime = torch.arange(nSample, dtype=x.real.dtype) / sampFreq
        xLimit = (0, max((nSample - 1) / sampFreq, 1 / sampFreq))
    freq = torch.linspace(0, fsPlot / 2, spec.shape[0], dtype=x.real.dtype)

    for index in range(nSignal):
        spectrum = 20 * torch.log10(torch.abs(torch.fft.rfft(xForFft[:, index])).clamp_min(eps))
        specDb = 20 * torch.log10(torch.abs(spec[:, :, index]).clamp_min(eps))
        vmax = float(specDb.max().item()) - trunc
        vmin = vmax - dynamicRange

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
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("isDraw=True requires matplotlib. Install it with `pip install matplotlib`.") from exc

    costCpu = cost.detach().cpu()
    # [bssAuxFdica.m:359] figure; plot(0:nIter, cost);
    fig, ax = plt.subplots()
    ax.plot(range(nIter + 1), costCpu[: nIter + 1].numpy())
    # [bssAuxFdica.m:361] xlabel("Number of iterations"); ylabel("Value of cost function");
    ax.set_xlabel("Number of iterations")
    ax.set_ylabel("Value of cost function")
    # [bssAuxFdica.m:362] grid on;
    ax.grid(True)
    ax.tick_params(labelsize=12)


def local_showPlots() -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("isDraw=True requires matplotlib. Install it with `pip install matplotlib`.") from exc

    plt.show(block=True)


# [bssAuxFdica.m:221] function [Y, W, cost] = local_auxFdica(X, nIter, srcModel, isDraw)
def local_auxFdica(
    X: torch.Tensor,
    nIter: int,
    srcModel: str,
    isDraw: bool = False,
    verbose: bool = False,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    # [bssAuxFdica.m:236] [I, J, M] = size(X, [1,2,3]); % nFreq x nTime x nCh
    # [bssAuxFdica.m:237] N = M; % number of sources
    nFreq, nFrame, nCh = X.shape
    # [bssAuxFdica.m:238] E = repmat(eye(M), [1, 1, I]);
    # [bssAuxFdica.m:239] W = E; % initial demixing matrix (N x M x I)
    W = torch.eye(nCh, dtype=X.dtype, device=X.device).repeat(nFreq, 1, 1)
    # [bssAuxFdica.m:240] Y = X; % initial estimated spectrogram
    Y = X.clone()
    # [bssAuxFdica.m:244] cost = zeros(nIter, 1);
    cost = torch.zeros(nIter + 1 if isDraw else nIter, dtype=X.real.dtype, device=X.device)
    if isDraw:
        # [bssAuxFdica.m:246] cost(1,1) = local_calcFdicaCost(Yp, W, srcModel, I, J);
        cost[0] = local_calcFdicaCost(Y, W, srcModel)
    # [bssAuxFdica.m:254] Rp = max(abs(Yp), 10000*eps);
    threshold = 10000 * torch.finfo(X.real.dtype).eps
    # [bssAuxFdica.m:238] E = repmat(eye(M), [1, 1, I]);
    eye = torch.eye(nCh, dtype=X.dtype, device=X.device)
    # [bssAuxFdica.m:251] for iIter = 1:nIter
    for iIter in range(nIter):
        iterationStart = time.perf_counter()
        if verbose:
            print(f"[FDICA] 学習 {iIter + 1}/{nIter} 開始", flush=True)
        # [bssAuxFdica.m:254] Rp = max(abs(Yp), 10000*eps);
        # [bssAuxFdica.m:255] elseif srcModel == "TVG"
        # [bssAuxFdica.m:256] Rp = max(abs(Yp).^2, 10000*eps);
        radius = torch.abs(Y) if srcModel == "LAP" else torch.abs(Y) ** 2
        # [bssAuxFdica.m:259] invRp = 1./Rp; % N x J x I
        invRadius = radius.clamp_min(threshold).reciprocal()
        # [bssAuxFdica.m:260] for n = 1:N
        for n in range(nCh):
            # [bssAuxFdica.m:262] Vk = pagemtimes(D.*Xp, Xph)/J; % M x M x I, pagewise matrix multiplication ((D(:,:,i).*Xp(:,:,i))*Xp(:,:,i)'/J)
            for f in range(nFreq):
                # [bssAuxFdica.m:241] Xp = permute(X, [3, 2, 1]); % M x J x I
                Xf = X[f].mT
                # [bssAuxFdica.m:261] D = repmat(invRp(n, :, :), [M, 1, 1]); % M x J x I
                # [bssAuxFdica.m:262] Vk = pagemtimes(D.*Xp, Xph)/J; % M x M x I, pagewise matrix multiplication ((D(:,:,i).*Xp(:,:,i))*Xp(:,:,i)'/J)
                Vn = (Xf * invRadius[f, :, n][None, :]) @ Xf.mH / nFrame
                meanPower = torch.real(torch.trace(Vn)) / nCh
                ridge = 100 * torch.finfo(X.real.dtype).eps * meanPower.clamp_min(1.0)
                VnReg = Vn + ridge * eye
                # [bssAuxFdica.m:263] wn = pagemldivide(pagemtimes(W, Vk), E(:, n, :)); % M x 1 x I, pagewise operation ((W(:,:,i)*Vk(:,:,i)) \ E(:, n, :))
                system = W[f] @ VnReg
                try:
                    wn = torch.linalg.solve(system, eye[:, n])
                except RuntimeError as exc:
                    if "singular" not in str(exc).lower():
                        raise
                    wn = torch.linalg.pinv(system) @ eye[:, n]
                # [bssAuxFdica.m:264] wn = wn ./ sqrt( pagemtimes(pagemtimes(wn, "ctranspose", Vk, "none"), wn) ); % M x 1 x I, pagewise operation (wn(:,:,i)/sqrt(wn(:,:,i)'*Vk(:,:,i)*wn(:,:,i)))
                norm = torch.sqrt(torch.real(wn.conj() @ VnReg @ wn).clamp_min(threshold))
                wn = wn / norm
                # [bssAuxFdica.m:267] W(n, :, :) = wnh;
                W[f, n] = wn.conj()
                # [bssAuxFdica.m:266] Yp(n, :, :) = pagemtimes(wnh, Xp); % 1 x J x I, pagewise matrix multiplication (wnh(:,:,i)*Xp(:,:,i))
                Y[f, :, n] = X[f] @ wn.conj()
        if isDraw:
            # [bssAuxFdica.m:282] cost(iIter+1, 1) = local_calcFdicaCost(Yp, W, srcModel, I, J);
            cost[iIter + 1] = local_calcFdicaCost(Y, W, srcModel)
        if verbose:
            iterationTime = time.perf_counter() - iterationStart
            print(f"[FDICA] 学習 {iIter + 1}/{nIter} 完了 ({iterationTime:.2f}秒)", flush=True)
    # [bssAuxFdica.m:285] Y = permute(Yp, [3, 2, 1]);
    return Y, W.permute(1, 2, 0).contiguous(), cost


# [bssAuxFdica.m:303] function [fixY, fixW] = local_projectionBack(Y, S, W)
def local_projectionBack(
    Y: torch.Tensor, S: torch.Tensor, W: torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor]:
    if S.ndim == 2:
        S = S[:, :, None]
    nFreq, _, nSrc = Y.shape
    nRef = S.shape[2]
    fixY = torch.empty((nFreq, Y.shape[1], nSrc, nRef), dtype=Y.dtype, device=Y.device)
    fixW = torch.empty((nSrc, W.shape[1], nFreq, nRef), dtype=W.dtype, device=W.device)
    for f in range(nFreq):
        # [bssAuxFdica.m:321] Yp = permute(Y, [3, 2, 1]); % N x J x I
        Yf = Y[f].mT
        # [bssAuxFdica.m:324] Yph = pagectranspose(Yp); % J x N x I, pagewise Hermitian transpose (Yp')
        # [bssAuxFdica.m:325] YpYph = pagemtimes(Yp, Yph); % N x N x I, pagewise matrix multiplication (Yp*Yp')
        # [bssAuxFdica.m:326] YphOnYpYph = pagemrdivide(Yph, YpYph); % J x N x I, pagewise matrix right-division (Yp'/(Yp*Yp'))
        # [bssAuxFdica.m:327] A = pagemtimes(Sp, YphOnYpYph); % 1 x N x I or M x N x I, pagewise matrix multiplication (Sp * Yp'/(Yp*Yp'))
        A = S[f].mT @ Yf.mH @ torch.linalg.pinv(Yf @ Yf.mH)
        # [bssAuxFdica.m:330] fixY = Ap .* Ypp; % M x N x J x I, using implicit expansion
        # [bssAuxFdica.m:331] fixY = permute(fixY, [4, 3, 2, 1]); % I x J x N x M
        fixY[f] = (A[:, :, None] * Yf[None, :, :]).permute(2, 1, 0)
        # [bssAuxFdica.m:332] fixW = Ap .* Wp; % M x N x N x I, using implicit expansion
        # [bssAuxFdica.m:333] fixW = permute(fixW, [2, 3, 4, 1]); % N x N x I x M
        fixW[:, :, f, :] = (A[:, :, None] * W[:, :, f][None, :, :]).permute(1, 2, 0)
    return fixY.squeeze(-1) if nRef == 1 else fixY, fixW.squeeze(-1) if nRef == 1 else fixW


# [bssAuxFdica.m:155] if isFilt
def _time_domain_filter(obsSigInput: torch.Tensor, demixMat: torch.Tensor, fftSize: int) -> torch.Tensor:
    # [bssAuxFdica.m:158] W = cat(3, demixMatFix, flip(conj(demixMatFix(:, :, 2:end-1)), 3)); % produce beyond Nyquist components
    full = torch.cat((demixMat, torch.conj(torch.flip(demixMat[:, :, 1:-1], dims=(2,)))), dim=2)
    # [bssAuxFdica.m:159] demixFilt = real(ifft(W, fftSize, 3)); % fftSize x nSrc x nMic
    filt = torch.fft.ifft(full, n=fftSize, dim=2).real
    # [bssAuxFdica.m:160] demixFilt = circshift(demixFilt, fftSize/2+1, 3); % move peak to center by circular shifting
    filt = torch.roll(filt, shifts=fftSize // 2 + 1, dims=2)
    outputs = []
    # [bssAuxFdica.m:161] for iSrc = 1:nSrc
    for n in range(filt.shape[0]):
        total = None
        # [bssAuxFdica.m:162] for iCh = 1:nCh
        for ch in range(filt.shape[1]):
            # [bssAuxFdica.m:163] f = squeeze(demixFilt(iSrc, iCh, :));
            # [bssAuxFdica.m:164] tmp(:, iCh) = conv(obsSigInput(:, iCh), f); % linear convolution
            value = torch_functional.conv1d(
                obsSigInput[:, ch][None, None, :],
                torch.flip(filt[n, ch], dims=(0,))[None, None, :],
                padding=fftSize - 1,
            )[0, 0]
            # [bssAuxFdica.m:166] estSig(:, iSrc) = sum(tmp, 2);
            total = value if total is None else total + value
        outputs.append(total)
    # [bssAuxFdica.m:168] estSig(1:fftSize/2+1,:) = []; % cut initial components caused by group delay (circular shifting)
    return torch.stack(outputs, dim=1)[fftSize // 2 + 1 :]


# [bssAuxFdica.m:1] function [estSig, cost] = bssAuxFdica(obsSig, nSrc, args)
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
    # [bssAuxFdica.m:85] obsSig (:,:) double
    # [bssAuxFdica.m:86] nSrc (1,1) double {mustBeInteger, mustBePositive}
    # [bssAuxFdica.m:87] args.fftSize (1,1) double {mustBeInteger, mustBePositive} = 1024
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
    refs = [refMic] if isinstance(refMic, int) else list(refMic)
    if not refs or any(not isinstance(r, int) or r < 1 or r > obsSig.shape[1] for r in refs):
        raise ValueError("refMic contains an invalid one-based microphone index")
    if len(refs) != 1:
        raise ValueError("Python time-domain output currently requires one reference microphone")

    # [main.m:25] rng(seed);
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
    # [bssAuxFdica.m:121] F = DGTtool("windowName", "b", "windowLength", fftSize, "windowShift", shiftSize); % create DGTtool instance
    # [bssAuxFdica.m:122] obsSpec = F.DGT(obsSig); % STFT
    obsSpec = dgt_stft(obsSig, fftSize, shiftSize)
    if verbose:
        print(f"[FDICA] STFT完了: spectrum_shape={tuple(obsSpec.shape)}", flush=True)
        print("[FDICA] 白色化開始" if isWhiten else "[FDICA] 白色化をスキップ", flush=True)
    # [bssAuxFdica.m:126] if isWhiten
    # [bssAuxFdica.m:127] obsSpecInput = local_whitening(obsSpec, nSrc);
    # [bssAuxFdica.m:128] else
    # [bssAuxFdica.m:129] obsSpecInput = obsSpec(:, :, 1:nSrc); % discard unnecessary channels
    # [bssAuxFdica.m:130] end
    obsSpecInput = local_whitening(obsSpec, nSrc) if isWhiten else obsSpec[:, :, :nSrc]
    if verbose and isWhiten:
        print("[FDICA] 白色化完了", flush=True)
    if verbose:
        print(f"[FDICA] 学習開始: model={srcModel}, iterations={nIter}", flush=True)
    # [bssAuxFdica.m:133] [estSpecFdica, demixMat, cost] = local_auxFdica(obsSpecInput, nIter, srcModel, isDraw);
    estSpecFdica, demixMat, cost = local_auxFdica(obsSpecInput, nIter, srcModel, isDraw, verbose)
    if verbose:
        print("[FDICA] 学習完了", flush=True)
        print("[FDICA] Projection back開始", flush=True)
    # [bssAuxFdica.m:136] [estSpecFdicaFix, demixMatFix] = local_projectionBack(estSpecFdica, obsSpec(:,:,refMic), demixMat);
    fixed, demixFixed = local_projectionBack(estSpecFdica, obsSpec[:, :, refs[0] - 1], demixMat)
    if verbose:
        print("[FDICA] Projection back完了", flush=True)
        print(f"[FDICA] 置換問題の解決開始: solver={permSolver}", flush=True)

    # [bssAuxFdica.m:139] if permSolver == "none"
    if permSolver == "none":
        # [bssAuxFdica.m:140] estSpec = estSpecFdicaFix;
        estSpec = fixed
        # [bssAuxFdica.m:141] estPerm = repmat(1:nSrc, [nFreq, 1]);
        estPerm = torch.arange(nSrc, device=obsSig.device).repeat(fixed.shape[0], 1)
    elif permSolver == "COR":
        # [bssAuxFdica.m:143] [estSpec, estPerm] = permSolverCor(estSpecFdicaFix, args.isPowRatio, args.typeCor, args.deltaFreq, args.ratioFreq);
        estSpec, estPerm = permSolverCor(fixed, isPowRatio, typeCor, deltaFreq, ratioFreq)
    elif permSolver == "DOA":
        if micPos is None:
            raise ValueError("micPos is required for the DOA permutation solver")
        # [bssAuxFdica.m:145] [estSpec, estPerm] = permSolverDoa(demixMatFix, estSpecFdicaFix, args.micPos, args.sampFreq);
        estSpec, estPerm = permSolverDoa(demixFixed, fixed, torch.as_tensor(micPos), sampFreq, seed=seed)
    else:
        if srcSig is None or srcSig.ndim != 3:
            raise ValueError("srcSig=(sample, channel, source) is required for IPS")
        # [bssAuxFdica.m:147] srcSpect = F.DGT(squeeze(args.srcSig(:, args.refMic, :)));
        srcSpec = dgt_stft(srcSig[:, refs[0] - 1, :], fftSize, shiftSize)
        # [bssAuxFdica.m:148] [estSpec, estPerm] = permSolverIps(estSpecFdicaFix, srcSpect);
        estSpec, estPerm = permSolverIps(fixed, srcSpec)
    if verbose:
        print("[FDICA] 置換問題の解決完了", flush=True)
    # [bssAuxFdica.m:150] for iFreq = 1:nFreq
    # [bssAuxFdica.m:151] demixMatFix(:, :, iFreq) = demixMatFix(estPerm(iFreq, :), :, iFreq);
    # [bssAuxFdica.m:152] end
    demixFixed = torch.stack([demixFixed[estPerm[f], :, f] for f in range(fixed.shape[0])], dim=2)

    if verbose:
        method = "時間領域分離フィルタ" if isFilt else "逆STFT"
        print(f"[FDICA] 時間領域信号の生成開始: method={method}", flush=True)
    # [bssAuxFdica.m:155] if isFilt
    if isFilt:
        # [bssAuxFdica.m:157] obsSigInput = F.pinv(obsSpecInput); % observed signal
        obsSigInput = dgt_istft(obsSpecInput, fftSize, shiftSize)
        # [bssAuxFdica.m:173] estSig = estSig(1:sigLen, :);
        estSig = _time_domain_filter(obsSigInput, demixFixed, fftSize)[:sigLen]
    else:
        # [bssAuxFdica.m:171] estSig = F.pinv(estSpec);
        estSig = dgt_istft(estSpec, fftSize, shiftSize, length=sigLen)
    if verbose:
        print(f"[FDICA] 時間領域信号の生成完了: output_shape={tuple(estSig.shape)}", flush=True)
    if isDraw:
        if verbose:
            print("[FDICA] 描画開始: ウィンドウを閉じるまで処理を待機します", flush=True)
        # [bssAuxFdica.m:177] F.plot(obsSig, args.sampFreq); % observed signal
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
