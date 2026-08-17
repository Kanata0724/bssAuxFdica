"""周波数置換問題ソルバーの基本動作を確認するテスト。"""

import torch

from python_fdica.permSolverIps import permSolverIps


def test_ideal_permutation_solver() -> None:
    """1周波数だけ入れ替えた音源順序をIPSが正解へ戻せることを確認する。"""
    # 4周波数・6時間フレーム・2音源の正解スペクトログラムを作る。
    src = torch.randn(4, 6, 2, dtype=torch.complex128)
    mix = src.clone()
    # 周波数インデックス1だけ、音源0と音源1を意図的に交換する。
    mix[1] = src[1, :, [1, 0]]
    est, perm = permSolverIps(mix, src)
    # IPS後の推定値が正解と一致し、全周波数分の置換番号が返ることを確認する。
    torch.testing.assert_close(est, src)
    assert perm.shape == (4, 2)
