# 対応注: 学習用コピー：元ファイル python_fdica/permSolverCor.py（learning / b934fef 時点）。
"""相関を使って周波数ごとの音源順序をそろえるCOR置換問題ソルバー。

FDICAでは各周波数を別々に分離するため、ある周波数では1列目が話者A、
別の周波数では1列目が話者Bになる場合がある。本モジュールは時間方向の
パワー変化が同じ音源では似ることを利用し、全周波数で順序を統一する。
MATLAB版 ``permSolverCor.m`` に対応している。
"""

from __future__ import annotations

import itertools

import torch


# [permSolverCor.m:69] rhoGl = corr(vf, c); % correlation between feature and centroid vectors
# 対応注: MATLAB組込みcorrの役割をPythonの補助関数で実装。分母のeps下限はMATLAB呼び出しにはなく、定数列の扱いが異なる。
def _corr(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    """2つの行列の列どうしについてPearson相関係数を計算する。

    行は時間フレーム、列は音源を表す。戻り値の ``(i, j)`` は、
    ``a`` の音源iと ``b`` の音源jがどれだけ似ているかを表す。
    """
    # 平均を引き、各列を平均0の変動成分へ変換する。
    ac = a - a.mean(dim=0, keepdim=True)
    bc = b - b.mean(dim=0, keepdim=True)
    # 内積を各列の大きさで割ると、-1から1の相関係数になる。
    denom = torch.linalg.vector_norm(ac, dim=0)[:, None] * torch.linalg.vector_norm(bc, dim=0)[None, :]
    return (ac.mT @ bc) / denom.clamp_min(torch.finfo(a.dtype).eps)


def _local_frequency_set(f: int, nFreq: int, deltaFreq: int, ratioFreq: int) -> list[int]:
    """周波数 ``f`` の並びを決める際に比較する近傍・倍音周波数を返す。"""
    # まずfの前後deltaFreq個を近傍として集める。
    # [permSolverCor.m:145] adjSet = [f-delta:f-1, f+1:f+delta]; % set of adjacent local frequency, i.e., A(f)
    local = set(range(f - deltaFreq, f)) | set(range(f + 1, f + deltaFreq + 1))
    # [permSolverCor.m:147] for iRatio = 2:ratio % set of harmonic local frequency, i.e., H(f)
    for ratio in range(2, ratioFreq + 1):
        # 基本周波数の整数分の1と整数倍の周辺も、音色が共通する候補として使う。
        # MATLABのroundは0.5を絶対値の大きい側へ丸めるため、それも再現する。
        # [permSolverCor.m:148] harSet = [harSet, round(f/iRatio)-1:round(f/iRatio)+1, f*iRatio-1:f*iRatio+1];
        # 対応注: この行からlocal.updateまでで倍音近傍を構成。MATLABの1始まりへ変換してround相当を計算。
        lower = int(torch.floor(torch.tensor((f + 1) / ratio) + 0.5).item()) - 1
        harmonic = (f + 1) * ratio - 1
        local.update(range(lower - 1, lower + 2))
        local.update(range(harmonic - 1, harmonic + 2))
    # [permSolverCor.m:150] localSet = unique([adjSet, harSet]); % Union of A(f) and H(f) (and sorting)
    # [permSolverCor.m:151] localSet = localSet(localSet>=1 & localSet<=F); % frequency index must be in the range [1:nFreq]
    # 対応注: Pythonはf自身を除外するがMATLAB原文は明示的には除外しない。
    return sorted(i for i in local if 0 <= i < nFreq and i != f)


# [permSolverCor.m:1] function [est, perm] = permSolverCor(mix, isPowRatio, type, deltaFreq, ratioFreq)
def permSolverCor(
    mix: torch.Tensor,
    isPowRatio: bool = True,
    typeCor: str = "Gl+Lo",
    deltaFreq: int = 3,
    ratioFreq: int = 2,
) -> tuple[torch.Tensor, torch.Tensor]:
    """大域・局所相関を使い、各周波数の推定音源を同じ順番へそろえる。

    引数:
        mix: ``(frequency, frame, source)`` 形状の複素スペクトログラム。
        isPowRatio: Trueなら各音源のパワーを全音源パワーの合計で割る。
        typeCor: ``"Gl"`` は全周波数の大域相関、``"Lo"`` は近傍周波数の
            局所相関、``"Gl+Lo"`` は両方を使用する。
        deltaFreq: 局所相関で参照する前後の周波数ビン数。
        ratioFreq: 倍音関係として参照する最大倍率。

    戻り値:
        順序をそろえたスペクトログラムと、各周波数で使用した0始まりの
        置換番号 ``(frequency, source)``。
    """
    # [permSolverCor.m:31] mix (:,:,:) double
    # [permSolverCor.m:32] isPowRatio (1,1) logical = true
    # [permSolverCor.m:33] type {mustBeMember(type,{'Gl', 'Lo', 'Gl+Lo'})} = "Gl+Lo"
    # [permSolverCor.m:34] deltaFreq (1,1) {mustBeNonnegative} = 3
    # [permSolverCor.m:35] ratioFreq (1,1) {mustBeNonnegative} = 2
    # 対応注: 以下の入力検査と既定値はargumentsブロックに対応。ただし検査条件には差がある。
    if mix.ndim != 3 or not mix.is_complex():
        raise TypeError("mix must be a complex tensor shaped (frequency, frame, source)")
    if typeCor not in {"Gl", "Lo", "Gl+Lo"}:
        raise ValueError('typeCor must be "Gl", "Lo", or "Gl+Lo"')
    if not isinstance(deltaFreq, int) or deltaFreq < 0 or not isinstance(ratioFreq, int) or ratioFreq < 0:
        raise ValueError("deltaFreq and ratioFreq must be nonnegative integers")
    # [permSolverCor.m:37] [nFreq, nTime, nSrc] = size(mix, [1, 2, 3]);
    nFreq, _, nSrc = mix.shape
    # 複素振幅の絶対値二乗がパワーである。
    # [permSolverCor.m:42] v = abs(mix).^2 ./ sum(abs(mix).^2, 3); % power ratio Eq. (14), nFreq x nTime x nSrc
    # 対応注: 次行と合わせて特徴量を計算。Pythonは分母にeps下限を加える。
    power = torch.abs(mix) ** 2
    v = power / power.sum(dim=2, keepdim=True).clamp_min(torch.finfo(power.dtype).eps) if isPowRatio else torch.abs(mix)

    # 考えられる全音源順序を列挙する。MATLABのpermsと同じ探索順にするため
    # Pythonのitertools.permutationsの結果を逆順にする。
    # [permSolverCor.m:49] allPerm = perms((1:nSrc)); % all permutation patterns in nSrc sources case, nPerm x nSrc
    # 対応注: MATLABの候補順に合わせて逆順に列挙。番号は0始まり。
    allPerm = list(reversed(list(itertools.permutations(range(nSrc)))))
    # [permSolverCor.m:50] perm = repmat(allPerm(end,:), [nFreq, 1]); % initial permutation so that vPerm = v, nFreq x nSrc
    # 対応注: 次行のrepeatと合わせて初期置換を用意。
    identity = torch.arange(nSrc, device=mix.device)
    perm = identity.repeat(nFreq, 1)
    # [permSolverCor.m:51] vPerm = v; % initial permutation-fixed feature vector
    vPerm = v.clone()
    # [permSolverCor.m:54] while(true)
    while True:
        # 1周前と同じ置換になったら収束したと判断する。
        # [permSolverCor.m:58] permOld = perm; % permutation of previous iteration
        old = perm.clone()
        # centroidは、現在そろえられた全周波数の平均パターンである。
        # [permSolverCor.m:62] c = squeeze(mean(vPerm, 1)); % Eq. (17), nTime x nSrc
        # 対応注: MATLABのGlとGl+Lo分岐を共通化。
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
                # 対応注: 機能上の対応であり同値ではない。MATLABは周波数番号fを格納添字にして配列を拡張し得る。Pythonは候補だけをstackして平均する。
                rhoLo = torch.stack([_corr(vf, vPerm[g]) for g in localSet]).mean(dim=0)
            else:
                # 対応注: Python独自：候補周波数が空の場合は局所相関をゼロとする。
                rhoLo = torch.zeros((nSrc, nSrc), device=mix.device, dtype=v.dtype)
            # [permSolverCor.m:52] sumRho = zeros(nPerm, 1); % variable for storing cost in Eq. (18)
            costs = []
            rows = torch.arange(nSrc, device=mix.device)
            # [permSolverCor.m:72] for iPerm = 1:nPerm % calc Eq. (18) for all permutation patterns
            for p in allPerm:
                # 候補pで音源を対応付けたときの相関係数の合計を評価値にする。
                cols = torch.tensor(p, device=mix.device)
                value = torch.zeros((), device=mix.device, dtype=v.dtype)
                if rhoGl is not None:
                    # [permSolverCor.m:73] sumRho(iPerm, 1) = sum(diag(rhoGl(:, allPerm(iPerm, :)))); % diagonal elements of "rho(:, allPerm(iPerm, :))" are permuted combination
                    # 対応注: Pythonは対応する行・列の要素を選択して総和。
                    value = value + rhoGl[rows, cols].sum()
                if typeCor in {"Lo", "Gl+Lo"}:
                    # [permSolverCor.m:94] sumRho(iPerm, 1) = sum(diag(rhoLo));
                    # 対応注: Gl+Loでは大域相関の値へ加算。
                    value = value + rhoLo[rows, cols].sum()
                costs.append(value)
            # 最も相関が高くなる音源順序を、この周波数で採用する。
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
    # 求めた置換を元の複素スペクトログラムへ適用する。
    # [permSolverCor.m:137] est(iFreq, :, :) = mix(iFreq, :, perm(iFreq,:));
    est = torch.stack([mix[f, :, perm[f]] for f in range(nFreq)])
    return est, perm
