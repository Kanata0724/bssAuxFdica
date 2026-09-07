# 対応注: 学習用コピー：元ファイル python_fdica/tests/test_stft.py（learning / b934fef 時点）。
# 対応注: このテストファイルに直接対応するMATLABテストはない。Python独自の検証コードとして保持。
"""STFTと逆STFTが波形を正しく往復できることを確認するテスト。"""

import pytest
import torch

from python_fdica.stft import dgt_istft, dgt_stft


@pytest.mark.parametrize("shift_size", [8, 16])
def test_stft_istft_reconstructs_and_preserves_length(shift_size: int) -> None:
    """異なるフレームシフトでも、逆変換後の長さと値が元信号に一致する。"""
    # 乱数生成器のseedを固定し、毎回同じテスト信号を作る。
    generator = torch.Generator().manual_seed(4)
    x = torch.randn(257, 2, generator=generator, dtype=torch.float64)
    X = dgt_stft(x, fftSize=32, shiftSize=shift_size)
    y = dgt_istft(X, fftSize=32, shiftSize=shift_size, length=x.shape[0])
    # まず軸と長さを確認し、次に丸め誤差1e-10以内で全サンプルを比較する。
    assert y.shape == x.shape
    torch.testing.assert_close(y, x, atol=1e-10, rtol=1e-10)


def test_invalid_stft_arguments() -> None:
    """不正なFFT設定や整数波形を、分かりやすい例外として拒否する。"""
    # shiftSizeはfftSize以下でなければ、窓の間に未処理区間ができる。
    with pytest.raises(ValueError):
        dgt_stft(torch.ones(32, 2), fftSize=16, shiftSize=17)
    # STFTは実数の浮動小数点波形を入力として要求する。
    with pytest.raises(TypeError):
        dgt_stft(torch.ones(32, 2, dtype=torch.int64), fftSize=16, shiftSize=8)
