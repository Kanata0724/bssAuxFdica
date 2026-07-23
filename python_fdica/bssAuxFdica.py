"""Auxiliary-function FDICA corresponding to the MATLAB bssAuxFdica.m."""

from __future__ import annotations

from collections.abc import Sequence

import torch
import torch.nn.functional as torch_functional

from .permSolverCor import permSolverCor
from .permSolverDoa import permSolverDoa
from .permSolverIps import permSolverIps
from .stft import dgt_istft, dgt_stft


def local_whitening(X: torch.Tensor, nSrc: int) -> torch.Tensor:
    """Frequency-wise PCA whitening; X/Y are (frequency, frame, channel)."""
    nFreq, nFrame, _ = X.shape
    Y = torch.empty((nFreq, nFrame, nSrc), dtype=X.dtype, device=X.device)
    eps = torch.finfo(X.real.dtype).eps
    for iFreq in range(nFreq):
        Xi = X[iFreq].mT  # channel,frame
        covariance = Xi @ Xi.mH / nFrame
        eigenvalues, eigenvectors = torch.linalg.eigh(covariance)
        idx = torch.argsort(eigenvalues, descending=True)[:nSrc]
        dP = eigenvectors[:, idx]
        Y[iFreq] = ((dP.mH @ Xi) / torch.sqrt(eigenvalues[idx].clamp_min(eps))[:, None]).mT
    return Y


def local_calcFdicaCost(Y: torch.Tensor, W: torch.Tensor, srcModel: str) -> torch.Tensor:
    """Calculate the MATLAB objective from Y=(frequency,frame,source)."""
    nFrame = Y.shape[1]
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
    """Draw DGTtool-like spectrogram, spectrum, and waveform figures."""
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
    """Draw FDICA objective values with MATLAB-compatible axes."""
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
    """Show matplotlib figures until the user closes them."""
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
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Run LAP or TVG AuxFDICA on X=(frequency,frame,channel)."""
    nFreq, nFrame, nCh = X.shape
    W = torch.eye(nCh, dtype=X.dtype, device=X.device).repeat(nFreq, 1, 1)
    Y = X.clone()
    # Preserve MATLAB behavior: cost is only evaluated when plotting is requested.
    cost = torch.zeros(nIter + 1 if isDraw else nIter, dtype=X.real.dtype, device=X.device)
    if isDraw:
        cost[0] = local_calcFdicaCost(Y, W, srcModel)
    threshold = 10000 * torch.finfo(X.real.dtype).eps
    eye = torch.eye(nCh, dtype=X.dtype, device=X.device)
    for iIter in range(nIter):
        radius = torch.abs(Y) if srcModel == "LAP" else torch.abs(Y) ** 2
        invRadius = radius.clamp_min(threshold).reciprocal()  # frequency,frame,source
        for n in range(nCh):
            for f in range(nFreq):
                Xf = X[f].mT  # channel,frame
                Vn = (Xf * invRadius[f, :, n][None, :]) @ Xf.mH / nFrame
                wn = torch.linalg.solve(W[f] @ Vn, eye[:, n])
                norm = torch.sqrt(torch.real(wn.conj() @ Vn @ wn).clamp_min(threshold))
                wn = wn / norm
                W[f, n] = wn.conj()
                Y[f, :, n] = X[f] @ wn.conj()
        if isDraw:
            cost[iIter + 1] = local_calcFdicaCost(Y, W, srcModel)
    return Y, W.permute(1, 2, 0).contiguous(), cost


def local_projectionBack(
    Y: torch.Tensor, S: torch.Tensor, W: torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor]:
    """Project separated spectra onto reference channel(s)."""
    if S.ndim == 2:
        S = S[:, :, None]
    nFreq, _, nSrc = Y.shape
    nRef = S.shape[2]
    fixY = torch.empty((nFreq, Y.shape[1], nSrc, nRef), dtype=Y.dtype, device=Y.device)
    fixW = torch.empty((nSrc, W.shape[1], nFreq, nRef), dtype=W.dtype, device=W.device)
    for f in range(nFreq):
        Yf = Y[f].mT  # source,frame
        # A = S * Y^H / (Y*Y^H), using pinv for MATLAB right-division robustness.
        A = S[f].mT @ Yf.mH @ torch.linalg.pinv(Yf @ Yf.mH)  # reference,source
        fixY[f] = (A[:, :, None] * Yf[None, :, :]).permute(2, 1, 0)
        fixW[:, :, f, :] = (A[:, :, None] * W[:, :, f][None, :, :]).permute(1, 2, 0)
    return fixY.squeeze(-1) if nRef == 1 else fixY, fixW.squeeze(-1) if nRef == 1 else fixW


def _time_domain_filter(obsSigInput: torch.Tensor, demixMat: torch.Tensor, fftSize: int) -> torch.Tensor:
    """Apply MATLAB's optional full-spectrum, circularly shifted demixing FIR."""
    full = torch.cat((demixMat, torch.conj(torch.flip(demixMat[:, :, 1:-1], dims=(2,)))), dim=2)
    filt = torch.fft.ifft(full, n=fftSize, dim=2).real
    filt = torch.roll(filt, shifts=fftSize // 2 + 1, dims=2)
    outputs = []
    for n in range(filt.shape[0]):
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
) -> tuple[torch.Tensor, torch.Tensor]:
    """Separate a determined/overdetermined mixture with MATLAB-compatible FDICA.

    Args:
        obsSig: Real tensor shaped ``(sample, channel)``.
        nSrc: Number of sources. Internal axes are frequency, frame, source.
        srcModel: ``"LAP"`` or ``"TVG"``; this variable selects the update.
        refMic: MATLAB-compatible one-based reference microphone index/indices.
        permSolver: ``"none"``, ``"COR"``, ``"DOA"``, or ``"IPS"``.

    Returns:
        ``(estSig, cost)`` where ``estSig`` is ``(sample, source)`` for a
        scalar reference microphone. Cost follows MATLAB's plotting behavior.
    """
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

    torch.manual_seed(seed)
    sigLen = obsSig.shape[0]
    obsSpec = dgt_stft(obsSig, fftSize, shiftSize)
    obsSpecInput = local_whitening(obsSpec, nSrc) if isWhiten else obsSpec[:, :, :nSrc]
    estSpecFdica, demixMat, cost = local_auxFdica(obsSpecInput, nIter, srcModel, isDraw)
    fixed, demixFixed = local_projectionBack(estSpecFdica, obsSpec[:, :, refs[0] - 1], demixMat)

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
    demixFixed = torch.stack([demixFixed[estPerm[f], :, f] for f in range(fixed.shape[0])], dim=2)

    if isFilt:
        obsSigInput = dgt_istft(obsSpecInput, fftSize, shiftSize)
        estSig = _time_domain_filter(obsSigInput, demixFixed, fftSize)[:sigLen]
    else:
        estSig = dgt_istft(estSpec, fftSize, shiftSize, length=sigLen)
    if isDraw:
        local_plotSpectrogram(obsSig, sampFreq, fftSize, shiftSize, title="Observed signal")
        local_plotSpectrogram(dgt_istft(obsSpecInput, fftSize, shiftSize, length=sigLen), sampFreq, fftSize, shiftSize, title="FDICA input signal")
        local_plotSpectrogram(dgt_istft(estSpecFdica, fftSize, shiftSize, length=sigLen), sampFreq, fftSize, shiftSize, title="Estimated signal before projection back")
        if permSolver != "none":
            local_plotSpectrogram(dgt_istft(fixed, fftSize, shiftSize, length=sigLen), sampFreq, fftSize, shiftSize, title="Estimated signal before permutation solver")
        local_plotSpectrogram(estSig, sampFreq, fftSize, shiftSize, title="Estimated signal")
        local_plotCost(cost, nIter)
        local_showPlots()
    return estSig, cost
