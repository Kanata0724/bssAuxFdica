"""MATLAB版DGTtoolと互換性を持つSTFT・逆STFT処理。

時間波形を周波数領域へ変換する :func:`dgt_stft` と、その逆変換を行う
:func:`dgt_istft` を提供する。計算にはPyTorchを使っているため、これらの
関数を含む処理はオートグラドによるバックプロパゲーションが可能である。
"""

from __future__ import annotations

import math

import torch


def blackman_window(fftSize: int, *, device: torch.device, dtype: torch.dtype) -> torch.Tensor:
    """MATLAB版DGTtoolの ``windowName="b"`` に対応するBlackman窓を作る。

    ``periodic=True`` は、窓の両端を周期信号として扱うFFT用の定義である。
    ``device`` と ``dtype`` を入力信号に合わせることで、CPU/GPU間や
    float32/float64間の不要な変換を防ぐ。
    """
    return torch.blackman_window(fftSize, periodic=True, device=device, dtype=dtype)


def _check_transform_args(fftSize: int, shiftSize: int) -> None:
    """FFTサイズとフレームシフトがSTFTで使用可能か検査する。"""
    if not isinstance(fftSize, int) or fftSize <= 0:
        raise ValueError("fftSize must be a positive integer")
    if not isinstance(shiftSize, int) or shiftSize <= 0:
        raise ValueError("shiftSize must be a positive integer")
    if shiftSize > fftSize:
        raise ValueError("shiftSize must be less than or equal to fftSize")


def dgt_stft(obsSig: torch.Tensor, fftSize: int = 1024, shiftSize: int = 512) -> torch.Tensor:
    """時間波形をMATLAB版DGTtool互換の短時間フーリエ変換へ変換する。

    引数:
        obsSig: 実数の時間波形。形状は ``(sample, channel)`` または
            1チャンネルを表す ``(sample,)``。
        fftSize: 1フレームに含めるサンプル数。周波数分解能を決める。
        shiftSize: 隣接するフレーム間を何サンプルずらすかを表す。

    戻り値:
        複素スペクトログラム。形状は ``(frequency, frame, channel)``。
        実数信号では負の周波数が正の周波数から復元できるため、周波数数は
        ``fftSize // 2 + 1`` となる。
    """
    _check_transform_args(fftSize, shiftSize)
    if not isinstance(obsSig, torch.Tensor):
        raise TypeError("obsSig must be a torch.Tensor")
    if obsSig.ndim == 1:
        # 1次元波形も、以降の処理では「1チャンネルの2次元波形」として扱う。
        obsSig = obsSig[:, None]
    if obsSig.ndim != 2 or obsSig.shape[0] == 0:
        raise ValueError("obsSig must have shape (sample, channel) and be non-empty")
    if obsSig.is_complex() or not obsSig.dtype.is_floating_point:
        raise TypeError("obsSig must be a real floating-point tensor")

    # STFTの全フレームが同じ長さになるよう、末尾を0で埋める。
    sigLen, _ = obsSig.shape
    paddedLen = max(fftSize, math.ceil(sigLen / shiftSize) * shiftSize)
    if paddedLen % shiftSize:
        paddedLen = math.ceil(paddedLen / shiftSize) * shiftSize
    x = torch.nn.functional.pad(obsSig, (0, 0, 0, paddedLen - sigLen))

    # DGTtoolは先頭より前にも窓が重なる。信号末尾を先頭へ付け、周期的な
    # 境界として再現する（overlapは隣の窓と重なるサンプル数）。
    overlap = fftSize - shiftSize
    xPeriodic = torch.cat((x[-overlap:], x), dim=0) if overlap else x
    win = blackman_window(fftSize, device=x.device, dtype=x.dtype)
    # torch.stftは入力を (channel, sample) として受け取るので転置する。
    X = torch.stft(
        xPeriodic.T,
        n_fft=fftSize,
        hop_length=shiftSize,
        win_length=fftSize,
        window=win,
        center=False,
        normalized=False,
        onesided=True,
        return_complex=True,
    )
    # FDICA内で統一している (frequency, frame, channel) の順へ並べ替える。
    return X.permute(1, 2, 0).contiguous()


def dgt_istft(
    spec: torch.Tensor,
    fftSize: int = 1024,
    shiftSize: int = 512,
    *,
    length: int | None = None,
) -> torch.Tensor:
    """:func:`dgt_stft` の複素スペクトログラムを時間波形へ戻す。

    引数:
        spec: ``(frequency, frame, channel)`` 形状の複素テンソル。
        fftSize: STFT時と同じFFTサイズ。
        shiftSize: STFT時と同じフレームシフト。
        length: 必要な出力サンプル数。省略時はDGTtoolの擬似逆変換と同様に
            ``frame * shiftSize`` サンプルを返す。

    戻り値:
        ``(sample, channel)`` 形状の実数時間波形。
    """
    _check_transform_args(fftSize, shiftSize)
    if not isinstance(spec, torch.Tensor) or spec.ndim != 3:
        raise ValueError("spec must have shape (frequency, frame, channel)")
    if not spec.is_complex():
        raise TypeError("spec must be complex-valued")
    if spec.shape[0] != fftSize // 2 + 1 or spec.shape[1] == 0:
        raise ValueError("spec has an incompatible frequency or frame dimension")
    # STFT時に0埋めした長さ。lengthを指定した場合は最後に元の長さへ切る。
    paddedLen = spec.shape[1] * shiftSize
    if length is not None and (not isinstance(length, int) or length <= 0 or length > paddedLen):
        raise ValueError("length must be a positive integer no greater than frame * shiftSize")

    # torch.istftの通常の線形OLA（重ね合わせ加算）をDGTtoolの周期的OLAへ
    # 合わせるため、先頭と末尾のフレームを循環コピーする。
    overlap = fftSize - shiftSize
    nWrap = max(1, math.ceil(overlap / shiftSize))
    frameIndex = torch.arange(-nWrap, spec.shape[1] + nWrap, device=spec.device) % spec.shape[1]
    wrapped = spec[:, frameIndex, :].permute(2, 0, 1).contiguous()
    real_dtype = spec.real.dtype
    win = blackman_window(fftSize, device=spec.device, dtype=real_dtype)
    # torch.istftが期待する軸順 (channel, frequency, frame) で逆変換する。
    raw = torch.istft(
        wrapped,
        n_fft=fftSize,
        hop_length=shiftSize,
        win_length=fftSize,
        window=win,
        center=True,
        normalized=False,
        onesided=True,
        length=(wrapped.shape[-1] - 1) * shiftSize,
    )
    # 循環コピーとcenter=Trueによって追加された前後部分を取り除く。
    offset = nWrap * shiftSize + overlap - fftSize // 2
    x = raw[:, offset : offset + paddedLen].T.contiguous()
    return x[:length] if length is not None else x
