from __future__ import annotations

import itertools

import torch


# [permSolverIps.m:1] function [est, perm] = permSolverIps(mix, src)
def permSolverIps(mix: torch.Tensor, src: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    # [permSolverIps.m:25] if ~isequal(size(mix), size(src)); error("Sizes of 'mix' and 'src' must be equal.\n"); end
    if mix.ndim != 3 or src.ndim != 3 or mix.shape != src.shape:
        raise ValueError("mix and src must have the same (frequency, frame, source) shape")
    if not mix.is_complex() or not src.is_complex():
        raise TypeError("mix and src must be complex-valued")
    # [permSolverIps.m:22] [nFreq, nTime, nSrc] = size(mix, [1, 2, 3]);
    nFreq, _, nSrc = mix.shape
    # [permSolverIps.m:30] permAll = perms(1:nSrc);
    allPerm = list(itertools.permutations(range(nSrc)))
    # [permSolverIps.m:28] est = zeros(nFreq, nTime, nSrc);
    # [permSolverIps.m:29] perm = zeros(nFreq, nSrc);
    est = torch.empty_like(mix)
    perm = torch.empty((nFreq, nSrc), dtype=torch.long, device=mix.device)
    # [permSolverIps.m:32] for iFreq = 1:nFreq
    for iFreq in range(nFreq):
        errors = []
        # [permSolverIps.m:33] for iPerm = 1:nPerm
        for p in allPerm:
            # [permSolverIps.m:37] err(iPerm) = err(iPerm) + sum(abs(mix(iFreq, :, mixInd) - src(iFreq, :, iSrc)).^2, "all");
            errors.append(torch.sum(torch.abs(mix[iFreq, :, list(p)] - src[iFreq]) ** 2))
        # [permSolverIps.m:40] [~, ind] = min(err);
        best = int(torch.argmin(torch.stack(errors)).item())
        # [permSolverIps.m:41] perm(iFreq, :) = permAll(ind, :);
        perm[iFreq] = torch.tensor(allPerm[best], device=mix.device)
        # [permSolverIps.m:42] est(iFreq, :, :) = mix(iFreq, :, perm(iFreq,:));
        est[iFreq] = mix[iFreq, :, perm[iFreq]]
    return est, perm
