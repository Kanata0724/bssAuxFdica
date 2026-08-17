"""正解音源を利用して周波数ごとの音源順序をそろえるIPS置換問題ソルバー。

FDICAは周波数ごとに独立して分離するため、周波数ビンごとに「音源1」と
「音源2」の並びが入れ替わることがある。IPSは正解スペクトログラムとの
誤差が最小になる並びを総当たりで選ぶ、評価・実験用のoracle手法である。
"""

from __future__ import annotations

import itertools

import torch


def permSolverIps(mix: torch.Tensor, src: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    """各周波数の推定音源を、正解音源 ``src`` と同じ順番へ並べ替える。

    ``mix`` と ``src`` はどちらも ``(frequency, frame, source)`` 形状である。
    戻り値は、並べ替え後のスペクトログラムと0始まりの置換番号である。
    """
    if mix.ndim != 3 or src.ndim != 3 or mix.shape != src.shape:
        raise ValueError("mix and src must have the same (frequency, frame, source) shape")
    if not mix.is_complex() or not src.is_complex():
        raise TypeError("mix and src must be complex-valued")
    nFreq, _, nSrc = mix.shape
    # 音源数が2なら (0, 1) と (1, 0) の2通りを作る。
    allPerm = list(itertools.permutations(range(nSrc)))
    est = torch.empty_like(mix)
    perm = torch.empty((nFreq, nSrc), dtype=torch.long, device=mix.device)
    for iFreq in range(nFreq):
        # 各周波数について、すべての並び順と正解との二乗誤差を比較する。
        errors = []
        for p in allPerm:
            errors.append(torch.sum(torch.abs(mix[iFreq, :, list(p)] - src[iFreq]) ** 2))
        # 誤差が最小の並びを、この周波数の音源順序として採用する。
        best = int(torch.argmin(torch.stack(errors)).item())
        perm[iFreq] = torch.tensor(allPerm[best], device=mix.device)
        est[iFreq] = mix[iFreq, :, perm[iFreq]]
    return est, perm
