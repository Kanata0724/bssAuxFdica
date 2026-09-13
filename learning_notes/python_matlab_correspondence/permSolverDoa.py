from __future__ import annotations

import math

import torch


# [permSolverDoa.m:1] function [est, perm] = permSolverDoa(demixMat, mix, micPos, sampFreq)
def permSolverDoa(
    demixMat: torch.Tensor,
    mix: torch.Tensor,
    micPos: torch.Tensor,
    sampFreq: float,
    *,
    seed: int = 0,
) -> tuple[torch.Tensor, torch.Tensor]:
    del seed
    # [permSolverDoa.m:29] arguments
    if mix.ndim != 3 or demixMat.ndim != 3:
        raise ValueError("mix and demixMat must be three-dimensional")
    if not mix.is_complex() or not demixMat.is_complex():
        raise TypeError("mix and demixMat must be complex-valued")
    # [permSolverDoa.m:35] [nFreq, nTime, nSrc] = size(mix, [1, 2, 3]);
    nFreq, _, nSrc = mix.shape
    if nSrc != 2 or demixMat.shape[:2] != (2, 2) or demixMat.shape[2] != nFreq:
        raise ValueError("DOA solver is implemented only for two sources/channels")
    micPos = torch.as_tensor(micPos, dtype=mix.real.dtype, device=mix.device).flatten()
    if micPos.numel() != 2 or torch.any(micPos < 0) or sampFreq <= 0:
        raise ValueError("micPos must contain two nonnegative positions and sampFreq must be positive")
    spacing = torch.abs(micPos[0] - micPos[1])
    if spacing == 0:
        raise ValueError("microphone positions must be distinct")

    # [permSolverDoa.m:49] mixMat(:, :, iFreq) = inv(demixMat(:,:,iFreq)); % mixMat(:, n, iFreq) is a steering vector of n-th source
    W = demixMat.permute(2, 0, 1)
    A = torch.linalg.inv(W)
    # [permSolverDoa.m:44] freqAx = linspace(0, sampFreq/2, nFreq).'; % nFreq x 1
    freqAx = torch.linspace(0, sampFreq / 2, nFreq, device=mix.device, dtype=mix.real.dtype)
    # [permSolverDoa.m:45] sinDoa = zeros(nFreq, nSrc);
    sinDoa = torch.full((nFreq, 2), float("nan"), device=mix.device, dtype=mix.real.dtype)
    denom = 2 * math.pi * freqAx[1:, None] * spacing
    # [permSolverDoa.m:56] srcAngle = permute(angle(mixMat(1, iSrc, nonZero)./mixMat(2, iSrc, nonZero)), [3, 1, 2]); % nFreq x nCh x nSrc
    phase = torch.angle(A[1:, 0, :] / A[1:, 1, :])
    # [permSolverDoa.m:57] sinDoa(nonZero, iSrc)  = srcAngle ./ (2*pi*freqAx(nonZero) * abs(micPos(1)-micPos(2))) * soundSpd;
    sinDoa[1:] = phase / denom * 340.0
    # [permSolverDoa.m:59] isValid(nonZero, iSrc) = abs(sinDoa(nonZero, iSrc)) < 1;
    valid = torch.isfinite(sinDoa) & (torch.abs(sinDoa) < 1)
    # [permSolverDoa.m:58] doa(nonZero, iSrc)     = real(asin(sinDoa(nonZero, iSrc))) / pi*180;
    doa = torch.rad2deg(torch.asin(torch.clamp(sinDoa, -1, 1)))
    # [permSolverDoa.m:64] validDOA = doa(isValid);
    values = doa[valid]
    if values.numel() < 2:
        raise ValueError("not enough valid DOA estimates for clustering")

    # [permSolverDoa.m:65] [~, centroids] = kmeans(validDOA, nSrc);
    centroids = torch.stack((values.min(), values.max()))
    for _ in range(100):
        labels = torch.argmin(torch.abs(values[:, None] - centroids[None, :]), dim=1)
        updated = torch.stack([values[labels == k].mean() if torch.any(labels == k) else centroids[k] for k in range(2)])
        if torch.allclose(updated, centroids):
            break
        centroids = updated
    # [permSolverDoa.m:66] boundary = mean(centroids); % boundary DOA of clusters
    boundary = centroids.mean()
    # [permSolverDoa.m:72] if doa(iFreq, 1) <= boundary && doa(iFreq, 2) >= boundary
    # [permSolverDoa.m:73] perm(iFreq, :) = [1, 2];
    # [permSolverDoa.m:74] else
    # [permSolverDoa.m:75] perm(iFreq, :) = [2, 1];
    # [permSolverDoa.m:76] end
    identity = (doa[:, 0] <= boundary) & (doa[:, 1] >= boundary)
    perm = torch.where(
        identity[:, None],
        torch.tensor([0, 1], device=mix.device),
        torch.tensor([1, 0], device=mix.device),
    )
    # [permSolverDoa.m:77] est(iFreq, :, :) = mix(iFreq, :, perm(iFreq, :));
    est = torch.stack([mix[f, :, perm[f]] for f in range(nFreq)])
    return est, perm
