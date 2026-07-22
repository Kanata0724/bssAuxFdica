"""DOA frequency permutation solver corresponding to permSolverDoa.m."""

from __future__ import annotations

import math

import torch


def permSolverDoa(
    demixMat: torch.Tensor,
    mix: torch.Tensor,
    micPos: torch.Tensor,
    sampFreq: float,
    *,
    seed: int = 0,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Align a two-source estimate by steering-vector direction of arrival.

    ``demixMat`` is ``(source, channel, frequency)`` and ``mix`` is
    ``(frequency, frame, source)``. ``seed`` is accepted for API-level
    reproducibility; the deterministic one-dimensional clustering needs no RNG.
    """
    del seed
    if mix.ndim != 3 or demixMat.ndim != 3:
        raise ValueError("mix and demixMat must be three-dimensional")
    if not mix.is_complex() or not demixMat.is_complex():
        raise TypeError("mix and demixMat must be complex-valued")
    nFreq, _, nSrc = mix.shape
    if nSrc != 2 or demixMat.shape[:2] != (2, 2) or demixMat.shape[2] != nFreq:
        raise ValueError("DOA solver is implemented only for two sources/channels")
    micPos = torch.as_tensor(micPos, dtype=mix.real.dtype, device=mix.device).flatten()
    if micPos.numel() != 2 or torch.any(micPos < 0) or sampFreq <= 0:
        raise ValueError("micPos must contain two nonnegative positions and sampFreq must be positive")
    spacing = torch.abs(micPos[0] - micPos[1])
    if spacing == 0:
        raise ValueError("microphone positions must be distinct")

    W = demixMat.permute(2, 0, 1)  # frequency,source,channel
    A = torch.linalg.inv(W)  # frequency,channel,source
    freqAx = torch.linspace(0, sampFreq / 2, nFreq, device=mix.device, dtype=mix.real.dtype)
    sinDoa = torch.full((nFreq, 2), float("nan"), device=mix.device, dtype=mix.real.dtype)
    denom = 2 * math.pi * freqAx[1:, None] * spacing
    phase = torch.angle(A[1:, 0, :] / A[1:, 1, :])
    sinDoa[1:] = phase / denom * 340.0
    valid = torch.isfinite(sinDoa) & (torch.abs(sinDoa) < 1)
    doa = torch.rad2deg(torch.asin(torch.clamp(sinDoa, -1, 1)))
    values = doa[valid]
    if values.numel() < 2:
        raise ValueError("not enough valid DOA estimates for clustering")

    # Deterministic 1-D k-means with extreme values as initial centroids.
    centroids = torch.stack((values.min(), values.max()))
    for _ in range(100):
        labels = torch.argmin(torch.abs(values[:, None] - centroids[None, :]), dim=1)
        updated = torch.stack([values[labels == k].mean() if torch.any(labels == k) else centroids[k] for k in range(2)])
        if torch.allclose(updated, centroids):
            break
        centroids = updated
    boundary = centroids.mean()
    identity = (doa[:, 0] <= boundary) & (doa[:, 1] >= boundary)
    perm = torch.where(
        identity[:, None],
        torch.tensor([0, 1], device=mix.device),
        torch.tensor([1, 0], device=mix.device),
    )
    est = torch.stack([mix[f, :, perm[f]] for f in range(nFreq)])
    return est, perm
