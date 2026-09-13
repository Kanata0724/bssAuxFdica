from __future__ import annotations

import itertools

import torch


# [permSolverCor.m:69] rhoGl = corr(vf, c); % correlation between feature and centroid vectors
def _corr(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    ac = a - a.mean(dim=0, keepdim=True)
    bc = b - b.mean(dim=0, keepdim=True)
    denom = torch.linalg.vector_norm(ac, dim=0)[:, None] * torch.linalg.vector_norm(bc, dim=0)[None, :]
    return (ac.mT @ bc) / denom.clamp_min(torch.finfo(a.dtype).eps)


def _local_frequency_set(f: int, nFreq: int, deltaFreq: int, ratioFreq: int) -> list[int]:
    # [permSolverCor.m:145] adjSet = [f-delta:f-1, f+1:f+delta]; % set of adjacent local frequency, i.e., A(f)
    local = set(range(f - deltaFreq, f)) | set(range(f + 1, f + deltaFreq + 1))
    # [permSolverCor.m:147] for iRatio = 2:ratio % set of harmonic local frequency, i.e., H(f)
    for ratio in range(2, ratioFreq + 1):
        # [permSolverCor.m:148] harSet = [harSet, round(f/iRatio)-1:round(f/iRatio)+1, f*iRatio-1:f*iRatio+1];
        lower = int(torch.floor(torch.tensor((f + 1) / ratio) + 0.5).item()) - 1
        harmonic = (f + 1) * ratio - 1
        local.update(range(lower - 1, lower + 2))
        local.update(range(harmonic - 1, harmonic + 2))
    # [permSolverCor.m:150] localSet = unique([adjSet, harSet]); % Union of A(f) and H(f) (and sorting)
    # [permSolverCor.m:151] localSet = localSet(localSet>=1 & localSet<=F); % frequency index must be in the range [1:nFreq]
    return sorted(i for i in local if 0 <= i < nFreq and i != f)


# [permSolverCor.m:1] function [est, perm] = permSolverCor(mix, isPowRatio, type, deltaFreq, ratioFreq)
def permSolverCor(
    mix: torch.Tensor,
    isPowRatio: bool = True,
    typeCor: str = "Gl+Lo",
    deltaFreq: int = 3,
    ratioFreq: int = 2,
) -> tuple[torch.Tensor, torch.Tensor]:
    # [permSolverCor.m:31] mix (:,:,:) double
    # [permSolverCor.m:32] isPowRatio (1,1) logical = true
    # [permSolverCor.m:33] type {mustBeMember(type,{'Gl', 'Lo', 'Gl+Lo'})} = "Gl+Lo"
    # [permSolverCor.m:34] deltaFreq (1,1) {mustBeNonnegative} = 3
    # [permSolverCor.m:35] ratioFreq (1,1) {mustBeNonnegative} = 2
    if mix.ndim != 3 or not mix.is_complex():
        raise TypeError("mix must be a complex tensor shaped (frequency, frame, source)")
    if typeCor not in {"Gl", "Lo", "Gl+Lo"}:
        raise ValueError('typeCor must be "Gl", "Lo", or "Gl+Lo"')
    if not isinstance(deltaFreq, int) or deltaFreq < 0 or not isinstance(ratioFreq, int) or ratioFreq < 0:
        raise ValueError("deltaFreq and ratioFreq must be nonnegative integers")
    # [permSolverCor.m:37] [nFreq, nTime, nSrc] = size(mix, [1, 2, 3]);
    nFreq, _, nSrc = mix.shape
    # [permSolverCor.m:42] v = abs(mix).^2 ./ sum(abs(mix).^2, 3); % power ratio Eq. (14), nFreq x nTime x nSrc
    power = torch.abs(mix) ** 2
    v = power / power.sum(dim=2, keepdim=True).clamp_min(torch.finfo(power.dtype).eps) if isPowRatio else torch.abs(mix)

    # [permSolverCor.m:49] allPerm = perms((1:nSrc)); % all permutation patterns in nSrc sources case, nPerm x nSrc
    allPerm = list(reversed(list(itertools.permutations(range(nSrc)))))
    # [permSolverCor.m:50] perm = repmat(allPerm(end,:), [nFreq, 1]); % initial permutation so that vPerm = v, nFreq x nSrc
    identity = torch.arange(nSrc, device=mix.device)
    perm = identity.repeat(nFreq, 1)
    # [permSolverCor.m:51] vPerm = v; % initial permutation-fixed feature vector
    vPerm = v.clone()
    # [permSolverCor.m:54] while(true)
    while True:
        # [permSolverCor.m:58] permOld = perm; % permutation of previous iteration
        old = perm.clone()
        # [permSolverCor.m:62] c = squeeze(mean(vPerm, 1)); % Eq. (17), nTime x nSrc
        centroid = vPerm.mean(dim=0) if typeCor in {"Gl", "Gl+Lo"} else None
        for f in range(nFreq):
            # [permSolverCor.m:66] vf = squeeze(v(iFreq, :, :)); % nTime x nSrc
            vf = v[f]
            # [permSolverCor.m:69] rhoGl = corr(vf, c); % correlation between feature and centroid vectors
            rhoGl = _corr(vf, centroid) if centroid is not None else None
            # [permSolverCor.m:86] localFreqSet = local_produceLocalFreqSet(iFreq, nFreq, deltaFreq, ratioFreq);
            localSet = _local_frequency_set(f, nFreq, deltaFreq, ratioFreq) if typeCor in {"Lo", "Gl+Lo"} else []
            if localSet:
                # [permSolverCor.m:158] rhoFreqwise(f, :, :) = corr(vf, vg); % correlation between feature and vg vectors
                rhoLo = torch.stack([_corr(vf, vPerm[g]) for g in localSet]).mean(dim=0)
            else:
                rhoLo = torch.zeros((nSrc, nSrc), device=mix.device, dtype=v.dtype)
            # [permSolverCor.m:52] sumRho = zeros(nPerm, 1); % variable for storing cost in Eq. (18)
            costs = []
            rows = torch.arange(nSrc, device=mix.device)
            # [permSolverCor.m:72] for iPerm = 1:nPerm % calc Eq. (18) for all permutation patterns
            for p in allPerm:
                cols = torch.tensor(p, device=mix.device)
                value = torch.zeros((), device=mix.device, dtype=v.dtype)
                if rhoGl is not None:
                    # [permSolverCor.m:73] sumRho(iPerm, 1) = sum(diag(rhoGl(:, allPerm(iPerm, :)))); % diagonal elements of "rho(:, allPerm(iPerm, :))" are permuted combination
                    value = value + rhoGl[rows, cols].sum()
                if typeCor in {"Lo", "Gl+Lo"}:
                    # [permSolverCor.m:94] sumRho(iPerm, 1) = sum(diag(rhoLo));
                    value = value + rhoLo[rows, cols].sum()
                costs.append(value)
            # [permSolverCor.m:163] [~, idx] = max(cost); % find index of maximum value
            best = int(torch.argmax(torch.stack(costs)).item())
            chosen = torch.tensor(allPerm[best], device=mix.device)
            # [permSolverCor.m:164] perm(f, :) = allPerm(idx, :); % permutation that maximizes Eq. (18)
            perm[f] = chosen
            # [permSolverCor.m:165] vPerm(f, :, :) = v(f, :, perm(f, :)); % update permutation-fixed v for calculating Eq. (17)
            vPerm[f] = v[f, :, chosen]
        # [permSolverCor.m:130] if all(permOld==perm, 'all'); break; end
        if torch.equal(old, perm):
            break
    # [permSolverCor.m:137] est(iFreq, :, :) = mix(iFreq, :, perm(iFreq,:));
    est = torch.stack([mix[f, :, perm[f]] for f in range(nFreq)])
    return est, perm
