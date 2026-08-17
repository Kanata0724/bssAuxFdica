"""到来方向（DOA）を使って周波数ごとの音源順序をそろえるソルバー。

2本のマイクに届く信号の位相差から、音がどちらの方向から来たかを推定する。
推定角度を2群へ分けることで、各周波数の2音源を同じ順番に並べる。
MATLAB版 ``permSolverDoa.m`` に対応している。
"""

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
    """2音源・2マイクの推定結果を、音源の到来方向に基づいてそろえる。

    ``demixMat`` は ``(source, channel, frequency)``、``mix`` は
    ``(frequency, frame, source)``。``micPos`` は2本のマイク位置、
    ``sampFreq`` はサンプリング周波数を表す。``seed`` は他ソルバーと
    APIをそろえるために受け取るが、この決定的な計算では乱数を使わない。
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

    # 分離行列Wの逆行列Aは、各音源が各マイクへ届く混合特性を表す。
    W = demixMat.permute(2, 0, 1)  # (frequency, source, channel)
    A = torch.linalg.inv(W)  # (frequency, channel, source)
    freqAx = torch.linspace(0, sampFreq / 2, nFreq, device=mix.device, dtype=mix.real.dtype)
    # 直流成分（0 Hz）は位相差から方向を求められないためNaNのままにする。
    sinDoa = torch.full((nFreq, 2), float("nan"), device=mix.device, dtype=mix.real.dtype)
    denom = 2 * math.pi * freqAx[1:, None] * spacing
    # 2マイク間の複素比の偏角が位相差になる。音速は340 m/sとしている。
    phase = torch.angle(A[1:, 0, :] / A[1:, 1, :])
    sinDoa[1:] = phase / denom * 340.0
    # arcsinへ渡せる範囲は-1から1なので、物理的に有効な推定だけを残す。
    valid = torch.isfinite(sinDoa) & (torch.abs(sinDoa) < 1)
    doa = torch.rad2deg(torch.asin(torch.clamp(sinDoa, -1, 1)))
    values = doa[valid]
    if values.numel() < 2:
        raise ValueError("not enough valid DOA estimates for clustering")

    # 角度を2音源に分ける1次元k-means。最小角と最大角を初期中心にするため、
    # 乱数に依存せず毎回同じ結果になる。
    centroids = torch.stack((values.min(), values.max()))
    for _ in range(100):
        labels = torch.argmin(torch.abs(values[:, None] - centroids[None, :]), dim=1)
        updated = torch.stack([values[labels == k].mean() if torch.any(labels == k) else centroids[k] for k in range(2)])
        if torch.allclose(updated, centroids):
            break
        centroids = updated
    # 2クラスタ中心の中間角を境界とし、左右どちらの音源かを決める。
    boundary = centroids.mean()
    identity = (doa[:, 0] <= boundary) & (doa[:, 1] >= boundary)
    perm = torch.where(
        identity[:, None],
        torch.tensor([0, 1], device=mix.device),
        torch.tensor([1, 0], device=mix.device),
    )
    # 周波数ごとに決定した順序をスペクトログラムへ適用する。
    est = torch.stack([mix[f, :, perm[f]] for f in range(nFreq)])
    return est, perm
