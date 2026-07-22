"""DGTtool-compatible STFT/ISTFT primitives implemented with PyTorch."""

from __future__ import annotations

import math

import torch


def blackman_window(fftSize: int, *, device: torch.device, dtype: torch.dtype) -> torch.Tensor:
    """Return the periodic Blackman window selected by DGTtool(windowName="b")."""
    return torch.blackman_window(fftSize, periodic=True, device=device, dtype=dtype)


def _check_transform_args(fftSize: int, shiftSize: int) -> None:
    if not isinstance(fftSize, int) or fftSize <= 0:
        raise ValueError("fftSize must be a positive integer")
    if not isinstance(shiftSize, int) or shiftSize <= 0:
        raise ValueError("shiftSize must be a positive integer")
    if shiftSize > fftSize:
        raise ValueError("shiftSize must be less than or equal to fftSize")


def dgt_stft(obsSig: torch.Tensor, fftSize: int = 1024, shiftSize: int = 512) -> torch.Tensor:
    """Compute DGTtool-compatible STFT.

    Args:
        obsSig: Real tensor shaped ``(sample, channel)`` or ``(sample,)``.

    Returns:
        Complex tensor shaped ``(frequency, frame, channel)``. The input is
        zero-padded at the end to a hop multiple and treated periodically at
        the analysis boundary, matching ``DGTtool.DGT``'s usual algorithm.
    """
    _check_transform_args(fftSize, shiftSize)
    if not isinstance(obsSig, torch.Tensor):
        raise TypeError("obsSig must be a torch.Tensor")
    if obsSig.ndim == 1:
        obsSig = obsSig[:, None]
    if obsSig.ndim != 2 or obsSig.shape[0] == 0:
        raise ValueError("obsSig must have shape (sample, channel) and be non-empty")
    if obsSig.is_complex() or not obsSig.dtype.is_floating_point:
        raise TypeError("obsSig must be a real floating-point tensor")

    sigLen, _ = obsSig.shape
    paddedLen = max(fftSize, math.ceil(sigLen / shiftSize) * shiftSize)
    if paddedLen % shiftSize:
        paddedLen = math.ceil(paddedLen / shiftSize) * shiftSize
    x = torch.nn.functional.pad(obsSig, (0, 0, 0, paddedLen - sigLen))

    # DGTtool frames start at -(windowLength-hop); prepend the periodic tail.
    overlap = fftSize - shiftSize
    xPeriodic = torch.cat((x[-overlap:], x), dim=0) if overlap else x
    win = blackman_window(fftSize, device=x.device, dtype=x.dtype)
    X = torch.stft(
        xPeriodic.T,
        n_fft=fftSize,
        hop_length=shiftSize,
        win_length=fftSize,
        window=win,
        center=False,
        normalized=False,
        onesided=True,
        return_complex=True,
    )
    return X.permute(1, 2, 0).contiguous()  # channel,freq,frame -> freq,frame,channel


def dgt_istft(
    spec: torch.Tensor,
    fftSize: int = 1024,
    shiftSize: int = 512,
    *,
    length: int | None = None,
) -> torch.Tensor:
    """Invert :func:`dgt_stft` using cyclic frame extension and ``torch.istft``.

    Args:
        spec: Complex tensor shaped ``(frequency, frame, channel)``.
        length: Optional requested sample count. Without it, returns
            ``frame * shiftSize`` samples, as DGTtool's pseudo-inverse does.

    Returns:
        Real tensor shaped ``(sample, channel)``.
    """
    _check_transform_args(fftSize, shiftSize)
    if not isinstance(spec, torch.Tensor) or spec.ndim != 3:
        raise ValueError("spec must have shape (frequency, frame, channel)")
    if not spec.is_complex():
        raise TypeError("spec must be complex-valued")
    if spec.shape[0] != fftSize // 2 + 1 or spec.shape[1] == 0:
        raise ValueError("spec has an incompatible frequency or frame dimension")
    paddedLen = spec.shape[1] * shiftSize
    if length is not None and (not isinstance(length, int) or length <= 0 or length > paddedLen):
        raise ValueError("length must be a positive integer no greater than frame * shiftSize")

    # Copies on both sides turn torch.istft's linear OLA into DGTtool's cyclic OLA.
    overlap = fftSize - shiftSize
    nWrap = max(1, math.ceil(overlap / shiftSize))
    frameIndex = torch.arange(-nWrap, spec.shape[1] + nWrap, device=spec.device) % spec.shape[1]
    wrapped = spec[:, frameIndex, :].permute(2, 0, 1).contiguous()
    real_dtype = spec.real.dtype
    win = blackman_window(fftSize, device=spec.device, dtype=real_dtype)
    raw = torch.istft(
        wrapped,
        n_fft=fftSize,
        hop_length=shiftSize,
        win_length=fftSize,
        window=win,
        center=True,
        normalized=False,
        onesided=True,
        length=(wrapped.shape[-1] - 1) * shiftSize,
    )
    offset = nWrap * shiftSize + overlap - fftSize // 2
    x = raw[:, offset : offset + paddedLen].T.contiguous()
    return x[:length] if length is not None else x
