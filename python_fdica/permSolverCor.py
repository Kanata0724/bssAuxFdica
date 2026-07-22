"""Correlation frequency permutation solver corresponding to permSolverCor.m."""

from __future__ import annotations

import itertools

import torch


def _corr(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    """MATLAB corr(a,b): columns are variables, rows are observations."""
    ac = a - a.mean(dim=0, keepdim=True)
    bc = b - b.mean(dim=0, keepdim=True)
    denom = torch.linalg.vector_norm(ac, dim=0)[:, None] * torch.linalg.vector_norm(bc, dim=0)[None, :]
    return (ac.mT @ bc) / denom.clamp_min(torch.finfo(a.dtype).eps)


def _local_frequency_set(f: int, nFreq: int, deltaFreq: int, ratioFreq: int) -> list[int]:
    local = set(range(f - deltaFreq, f)) | set(range(f + 1, f + deltaFreq + 1))
    for ratio in range(2, ratioFreq + 1):
        # MATLAB round is half-away-from-zero; indices here are nonnegative.
        lower = int(torch.floor(torch.tensor((f + 1) / ratio) + 0.5).item()) - 1
        harmonic = (f + 1) * ratio - 1
        local.update(range(lower - 1, lower + 2))
        local.update(range(harmonic - 1, harmonic + 2))
    return sorted(i for i in local if 0 <= i < nFreq and i != f)


def permSolverCor(
    mix: torch.Tensor,
    isPowRatio: bool = True,
    typeCor: str = "Gl+Lo",
    deltaFreq: int = 3,
    ratioFreq: int = 2,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Align frequency-wise source permutations using correlation clustering.

    Args:
        mix: Complex tensor ``(frequency, frame, source)``.
        typeCor: ``"Gl"``, ``"Lo"``, or ``"Gl+Lo"``.

    Returns:
        Aligned spectrogram and zero-based permutation array.
    """
    if mix.ndim != 3 or not mix.is_complex():
        raise TypeError("mix must be a complex tensor shaped (frequency, frame, source)")
    if typeCor not in {"Gl", "Lo", "Gl+Lo"}:
        raise ValueError('typeCor must be "Gl", "Lo", or "Gl+Lo"')
    if not isinstance(deltaFreq, int) or deltaFreq < 0 or not isinstance(ratioFreq, int) or ratioFreq < 0:
        raise ValueError("deltaFreq and ratioFreq must be nonnegative integers")
    nFreq, _, nSrc = mix.shape
    power = torch.abs(mix) ** 2
    v = power / power.sum(dim=2, keepdim=True).clamp_min(torch.finfo(power.dtype).eps) if isPowRatio else torch.abs(mix)

    # MATLAB perms is reverse lexicographic; its last row is the identity.
    allPerm = list(reversed(list(itertools.permutations(range(nSrc)))))
    identity = torch.arange(nSrc, device=mix.device)
    perm = identity.repeat(nFreq, 1)
    vPerm = v.clone()
    while True:
        old = perm.clone()
        centroid = vPerm.mean(dim=0) if typeCor in {"Gl", "Gl+Lo"} else None
        for f in range(nFreq):
            vf = v[f]
            rhoGl = _corr(vf, centroid) if centroid is not None else None
            localSet = _local_frequency_set(f, nFreq, deltaFreq, ratioFreq) if typeCor in {"Lo", "Gl+Lo"} else []
            if localSet:
                rhoLo = torch.stack([_corr(vf, vPerm[g]) for g in localSet]).mean(dim=0)
            else:
                rhoLo = torch.zeros((nSrc, nSrc), device=mix.device, dtype=v.dtype)
            costs = []
            rows = torch.arange(nSrc, device=mix.device)
            for p in allPerm:
                cols = torch.tensor(p, device=mix.device)
                value = torch.zeros((), device=mix.device, dtype=v.dtype)
                if rhoGl is not None:
                    value = value + rhoGl[rows, cols].sum()
                if typeCor in {"Lo", "Gl+Lo"}:
                    value = value + rhoLo[rows, cols].sum()
                costs.append(value)
            best = int(torch.argmax(torch.stack(costs)).item())
            chosen = torch.tensor(allPerm[best], device=mix.device)
            perm[f] = chosen
            vPerm[f] = v[f, :, chosen]
        if torch.equal(old, perm):
            break
    est = torch.stack([mix[f, :, perm[f]] for f in range(nFreq)])
    return est, perm
