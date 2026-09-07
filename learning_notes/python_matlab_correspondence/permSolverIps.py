# 対応注: 学習用コピー：元ファイル python_fdica/permSolverIps.py（learning / b934fef 時点）。
"""正解音源を利用して周波数ごとの音源順序をそろえるIPS置換問題ソルバー。

FDICAは周波数ごとに独立して分離するため、周波数ビンごとに「音源1」と
「音源2」の並びが入れ替わることがある。IPSは正解スペクトログラムとの
誤差が最小になる並びを総当たりで選ぶ、評価・実験用のoracle手法である。
"""

from __future__ import annotations

import itertools

import torch


# [permSolverIps.m:1] function [est, perm] = permSolverIps(mix, src)
def permSolverIps(mix: torch.Tensor, src: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    """各周波数の推定音源を、正解音源 ``src`` と同じ順番へ並べ替える。

    ``mix`` と ``src`` はどちらも ``(frequency, frame, source)`` 形状である。
    戻り値は、並べ替え後のスペクトログラムと0始まりの置換番号である。
    """
    # [permSolverIps.m:25] if ~isequal(size(mix), size(src)); error("Sizes of 'mix' and 'src' must be equal.\n"); end
    # 対応注: 次のis_complex検査もMATLABのisreal検査に対応。
    if mix.ndim != 3 or src.ndim != 3 or mix.shape != src.shape:
        raise ValueError("mix and src must have the same (frequency, frame, source) shape")
    if not mix.is_complex() or not src.is_complex():
        raise TypeError("mix and src must be complex-valued")
    # [permSolverIps.m:22] [nFreq, nTime, nSrc] = size(mix, [1, 2, 3]);
    nFreq, _, nSrc = mix.shape
    # 音源数が2なら (0, 1) と (1, 0) の2通りを作る。
    # [permSolverIps.m:30] permAll = perms(1:nSrc);
    # 対応注: PythonはMATLABと逆の候補順。誤差が同点の場合に選ばれる置換が異なり得る。
    allPerm = list(itertools.permutations(range(nSrc)))
    # [permSolverIps.m:28] est = zeros(nFreq, nTime, nSrc);
    # [permSolverIps.m:29] perm = zeros(nFreq, nSrc);
    # 対応注: Pythonは未初期化領域を確保し、後で全要素を代入。
    est = torch.empty_like(mix)
    perm = torch.empty((nFreq, nSrc), dtype=torch.long, device=mix.device)
    # [permSolverIps.m:32] for iFreq = 1:nFreq
    for iFreq in range(nFreq):
        # 各周波数について、すべての並び順と正解との二乗誤差を比較する。
        errors = []
        # [permSolverIps.m:33] for iPerm = 1:nPerm
        for p in allPerm:
            # [permSolverIps.m:37] err(iPerm) = err(iPerm) + sum(abs(mix(iFreq, :, mixInd) - src(iFreq, :, iSrc)).^2, "all");
            # 対応注: MATLABの音源ループと時間方向の総和をまとめて全要素の二乗誤差総和にする。
            errors.append(torch.sum(torch.abs(mix[iFreq, :, list(p)] - src[iFreq]) ** 2))
        # 誤差が最小の並びを、この周波数の音源順序として採用する。
        # [permSolverIps.m:40] [~, ind] = min(err);
        best = int(torch.argmin(torch.stack(errors)).item())
        # [permSolverIps.m:41] perm(iFreq, :) = permAll(ind, :);
        # 対応注: 番号はPythonが0始まり。
        perm[iFreq] = torch.tensor(allPerm[best], device=mix.device)
        # [permSolverIps.m:42] est(iFreq, :, :) = mix(iFreq, :, perm(iFreq,:));
        est[iFreq] = mix[iFreq, :, perm[iFreq]]
    return est, perm
