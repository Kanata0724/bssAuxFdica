"""Ideal (oracle) frequency permutation solver corresponding to permSolverIps.m."""

from __future__ import annotations

import itertools

import torch


def permSolverIps(mix: torch.Tensor, src: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    """Align each frequency bin to an oracle source spectrogram.

    Both inputs have shape ``(frequency, frame, source)``.
    """
    if mix.ndim != 3 or src.ndim != 3 or mix.shape != src.shape:
        raise ValueError("mix and src must have the same (frequency, frame, source) shape")
    if not mix.is_complex() or not src.is_complex():
        raise TypeError("mix and src must be complex-valued")
    nFreq, _, nSrc = mix.shape
    allPerm = list(itertools.permutations(range(nSrc)))
    est = torch.empty_like(mix)
    perm = torch.empty((nFreq, nSrc), dtype=torch.long, device=mix.device)
    for iFreq in range(nFreq):
        errors = []
        for p in allPerm:
            errors.append(torch.sum(torch.abs(mix[iFreq, :, list(p)] - src[iFreq]) ** 2))
        best = int(torch.argmin(torch.stack(errors)).item())
        perm[iFreq] = torch.tensor(allPerm[best], device=mix.device)
        est[iFreq] = mix[iFreq, :, perm[iFreq]]
    return est, perm
