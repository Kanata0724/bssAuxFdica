"""Python counterpart of main.m using repository dataset number 1."""

from __future__ import annotations

from pathlib import Path

import torch

from .bssAuxFdica import bssAuxFdica
from .example import read_pcm16_wav, write_pcm16_wav


def main() -> None:
    """Mix the two source images in dataset 1, separate them, and write WAVs."""
    repo = Path(__file__).resolve().parent.parent
    dataDir = repo / "dataset" / "dev1_female4_src_12_E2A_conv"
    names = [
        "dev1_female4_src_1_E2A_pos050130_mic2123_conv.wav",
        "dev1_female4_src_2_E2A_pos050130_mic2123_conv.wav",
    ]
    loaded = [read_pcm16_wav(dataDir / name) for name in names]
    images = [item[0] for item in loaded]
    sampleRate = loaded[0][1]
    if any(item[1] != sampleRate for item in loaded):
        raise ValueError("all source image files must have the same sample rate")
    obsSig = images[0] + images[1]
    peakVal = torch.max(torch.abs(obsSig))
    if peakVal > 1:
        # Intended main.m behavior (its clipping branch spells peakVal as peakValue).
        obsSig = 0.99 * obsSig / peakVal
        images = [0.99 * image / peakVal for image in images]
    estSig, _ = bssAuxFdica(
        obsSig,
        2,
        fftSize=4096,
        shiftSize=2048,
        nIter=50,
        isWhiten=False,
        srcModel="LAP",
        refMic=1,
        permSolver="COR",
        sampFreq=sampleRate,
        seed=1,
    )
    output = Path(__file__).resolve().parent / "output"
    output.mkdir(exist_ok=True)
    write_pcm16_wav(output / "data1_obs.wav", obsSig, sampleRate)
    write_pcm16_wav(output / "data1_src1.wav", images[0][:, :1], sampleRate)
    write_pcm16_wav(output / "data1_src2.wav", images[1][:, :1], sampleRate)
    write_pcm16_wav(output / "data1_est1.wav", estSig[:, :1], sampleRate)
    write_pcm16_wav(output / "data1_est2.wav", estSig[:, 1:2], sampleRate)
    print(f"wrote five WAV files under {output}; shape={tuple(estSig.shape)}")


if __name__ == "__main__":
    torch.set_default_dtype(torch.float64)
    main()
