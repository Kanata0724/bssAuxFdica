import pytest
import torch

from python_fdica import bssAuxFdica


def _mixture() -> torch.Tensor:
    t = torch.arange(320, dtype=torch.float64) / 8000
    src = torch.stack((torch.sin(2 * torch.pi * 440 * t), torch.sin(2 * torch.pi * 710 * t)), dim=1)
    return src @ torch.tensor([[1.0, 0.4], [0.3, 1.0]], dtype=torch.float64).mT


@pytest.mark.parametrize("model", ["LAP", "TVG"])
def test_lap_and_tvg_run_without_nonfinite_values(model: str) -> None:
    x = _mixture()
    y, cost = bssAuxFdica(
        x, 2, fftSize=32, shiftSize=16, nIter=2, isWhiten=False,
        srcModel=model, permSolver="none", seed=7,
    )
    assert y.shape == x.shape
    assert cost.shape == (2,)
    assert torch.isfinite(y).all()


def test_same_seed_is_reproducible() -> None:
    kwargs = dict(fftSize=32, shiftSize=16, nIter=2, isWhiten=False, srcModel="LAP", permSolver="none")
    first, _ = bssAuxFdica(_mixture(), 2, seed=11, **kwargs)
    second, _ = bssAuxFdica(_mixture(), 2, seed=11, **kwargs)
    torch.testing.assert_close(first, second, atol=0, rtol=0)


def test_time_domain_demixing_filter_runs() -> None:
    y, _ = bssAuxFdica(
        _mixture(), 2, fftSize=32, shiftSize=16, nIter=1,
        isWhiten=False, srcModel="LAP", permSolver="none", isFilt=True,
    )
    assert y.shape == _mixture().shape
    assert torch.isfinite(y).all()


@pytest.mark.parametrize(
    "kwargs",
    [
        {"nSrc": 3},
        {"nSrc": 2, "srcModel": "GAUSS"},
        {"nSrc": 2, "permSolver": "bad"},
        {"nSrc": 2, "nIter": 0},
        {"nSrc": 2, "refMic": 0},
    ],
)
def test_invalid_fdica_arguments(kwargs: dict) -> None:
    nSrc = kwargs.pop("nSrc")
    with pytest.raises(ValueError):
        bssAuxFdica(_mixture(), nSrc, fftSize=32, shiftSize=16, **kwargs)
