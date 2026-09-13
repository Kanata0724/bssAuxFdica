import pytest
import torch

from python_fdica.stft import dgt_istft, dgt_stft


@pytest.mark.parametrize("shift_size", [8, 16])
def test_stft_istft_reconstructs_and_preserves_length(shift_size: int) -> None:
    generator = torch.Generator().manual_seed(4)
    x = torch.randn(257, 2, generator=generator, dtype=torch.float64)
    X = dgt_stft(x, fftSize=32, shiftSize=shift_size)
    y = dgt_istft(X, fftSize=32, shiftSize=shift_size, length=x.shape[0])
    assert y.shape == x.shape
    torch.testing.assert_close(y, x, atol=1e-10, rtol=1e-10)


def test_invalid_stft_arguments() -> None:
    with pytest.raises(ValueError):
        dgt_stft(torch.ones(32, 2), fftSize=16, shiftSize=17)
    with pytest.raises(TypeError):
        dgt_stft(torch.ones(32, 2, dtype=torch.int64), fftSize=16, shiftSize=8)
