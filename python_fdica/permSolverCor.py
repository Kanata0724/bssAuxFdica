"""相関を使って周波数ごとの音源順序をそろえるCOR置換問題ソルバー。

FDICAでは各周波数を別々に分離するため、ある周波数では1列目が話者A、
別の周波数では1列目が話者Bになる場合がある。本モジュールは時間方向の
パワー変化が同じ音源では似ることを利用し、全周波数で順序を統一する。
MATLAB版 ``permSolverCor.m`` に対応している。
"""

from __future__ import annotations

import itertools

import torch


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
    local = set(range(f - deltaFreq, f)) | set(range(f + 1, f + deltaFreq + 1))
    for ratio in range(2, ratioFreq + 1):
        # 基本周波数の整数分の1と整数倍の周辺も、音色が共通する候補として使う。
        # MATLABのroundは0.5を絶対値の大きい側へ丸めるため、それも再現する。
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
    if mix.ndim != 3 or not mix.is_complex():
        raise TypeError("mix must be a complex tensor shaped (frequency, frame, source)")
    if typeCor not in {"Gl", "Lo", "Gl+Lo"}:
        raise ValueError('typeCor must be "Gl", "Lo", or "Gl+Lo"')
    if not isinstance(deltaFreq, int) or deltaFreq < 0 or not isinstance(ratioFreq, int) or ratioFreq < 0:
        raise ValueError("deltaFreq and ratioFreq must be nonnegative integers")
    nFreq, _, nSrc = mix.shape
    # 複素振幅の絶対値二乗がパワーである。
    power = torch.abs(mix) ** 2
    v = power / power.sum(dim=2, keepdim=True).clamp_min(torch.finfo(power.dtype).eps) if isPowRatio else torch.abs(mix)

    # 考えられる全音源順序を列挙する。MATLABのpermsと同じ探索順にするため
    # Pythonのitertools.permutationsの結果を逆順にする。
    allPerm = list(reversed(list(itertools.permutations(range(nSrc)))))
    identity = torch.arange(nSrc, device=mix.device)
    perm = identity.repeat(nFreq, 1)
    vPerm = v.clone()
    while True:
        # 1周前と同じ置換になったら収束したと判断する。
        old = perm.clone()
        # centroidは、現在そろえられた全周波数の平均パターンである。
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
                # 候補pで音源を対応付けたときの相関係数の合計を評価値にする。
                cols = torch.tensor(p, device=mix.device)
                value = torch.zeros((), device=mix.device, dtype=v.dtype)
                if rhoGl is not None:
                    value = value + rhoGl[rows, cols].sum()
                if typeCor in {"Lo", "Gl+Lo"}:
                    value = value + rhoLo[rows, cols].sum()
                costs.append(value)
            # 最も相関が高くなる音源順序を、この周波数で採用する。
            best = int(torch.argmax(torch.stack(costs)).item())
            chosen = torch.tensor(allPerm[best], device=mix.device)
            perm[f] = chosen
            vPerm[f] = v[f, :, chosen]
        if torch.equal(old, perm):
            break
    # 求めた置換を元の複素スペクトログラムへ適用する。
    est = torch.stack([mix[f, :, perm[f]] for f in range(nFreq)])
    return est, perm
