import torch

from python_fdica.permSolverIps import permSolverIps


def test_ideal_permutation_solver() -> None:
    src = torch.randn(4, 6, 2, dtype=torch.complex128)
    mix = src.clone()
    mix[1] = src[1, :, [1, 0]]
    est, perm = permSolverIps(mix, src)
    torch.testing.assert_close(est, src)
    assert perm.shape == (4, 2)
