"""Minimal tensor/WAV example for LAP and TVG FDICA."""

from __future__ import annotations

import argparse
import wave
from pathlib import Path

import torchaudio
import torch

from .bssAuxFdica import bssAuxFdica


def read_audio(path: str | Path, device: str = "cpu") -> tuple[torch.Tensor, int]:
    """Read an audio file with torchaudio as a ``(sample, channel)`` tensor."""
    values, sampleRate = torchaudio.load(str(path))
    return values.transpose(0, 1).to(device=device, dtype=torch.float64), sampleRate


def write_pcm16_wav(path: str | Path, signal: torch.Tensor, sampleRate: int) -> None:
    """Save ``(sample, channel)`` floating-point audio as a 16-bit PCM WAV."""
    signal16 = (signal.detach().cpu().clamp(-1, 1) * 32767).round().to(torch.int16).contiguous()
    with wave.open(str(path), "wb") as wav:
        wav.setnchannels(signal16.shape[1])
        wav.setsampwidth(2)
        wav.setframerate(sampleRate)
        wav.writeframes(bytes(signal16.view(torch.uint8).flatten().tolist()))


# Backward-compatible names used by main.py and existing callers.
read_pcm16_wav = read_audio


def make_example_mixture(sampleRate: int = 8000, seconds: float = 0.5) -> torch.Tensor:
    """Create a deterministic two-source/two-channel demonstration mixture."""
    t = torch.arange(int(sampleRate * seconds), dtype=torch.float64) / sampleRate
    sources = torch.stack((torch.sin(2 * torch.pi * 440 * t), torch.sin(2 * torch.pi * 733 * t)), dim=1)
    mixing = torch.tensor([[1.0, 0.45], [0.35, 1.0]], dtype=torch.float64)
    return sources @ mixing.mT


def main() -> None:
    """Load/create a mixture, run selected FDICA, and optionally save output."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, help="optional audio input supported by torchaudio")
    parser.add_argument("--output", type=Path, help="optional separated 16-bit PCM WAV")
    parser.add_argument("--model", choices=("LAP", "TVG"), default="LAP")
    parser.add_argument("--device", default="cpu", help='for example "cpu" or "cuda"')
    args = parser.parse_args()
    if args.input:
        obsSig, sampleRate = read_audio(args.input, args.device)
    else:
        sampleRate = 8000
        obsSig = make_example_mixture(sampleRate).to(args.device)
    estSig, _ = bssAuxFdica(
        obsSig,
        2,
        fftSize=256,
        shiftSize=128,
        nIter=10,
        srcModel=args.model,  # LAP/TVG is selected by this value.
        permSolver="COR",
        refMic=1,
        sampFreq=sampleRate,
        seed=1,
        isDraw=True
    )
    print(f"input={tuple(obsSig.shape)}, separated={tuple(estSig.shape)}, model={args.model}")
    if args.output:
        write_pcm16_wav(args.output, estSig, sampleRate)


if __name__ == "__main__":
    main()
