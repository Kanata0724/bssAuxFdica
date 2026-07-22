"""Minimal tensor/WAV example for LAP and TVG FDICA."""

from __future__ import annotations

import argparse
import wave
from pathlib import Path

import torch

from .bssAuxFdica import bssAuxFdica


def read_pcm16_wav(path: str | Path, device: str = "cpu") -> tuple[torch.Tensor, int]:
    """Read a 16-bit PCM WAV as a ``(sample, channel)`` float64 tensor."""
    with wave.open(str(path), "rb") as wav:
        if wav.getsampwidth() != 2:
            raise ValueError("the standard-library example supports 16-bit PCM WAV only")
        sampleRate = wav.getframerate()
        nChannel = wav.getnchannels()
        raw = wav.readframes(wav.getnframes())
    values = torch.frombuffer(bytearray(raw), dtype=torch.int16).reshape(-1, nChannel)
    return values.to(device=device, dtype=torch.float64) / 32768.0, sampleRate


def write_pcm16_wav(path: str | Path, signal: torch.Tensor, sampleRate: int) -> None:
    """Save ``(sample, channel)`` floating-point audio as 16-bit PCM WAV."""
    signal16 = (signal.detach().cpu().clamp(-1, 1) * 32767).round().to(torch.int16).contiguous()
    with wave.open(str(path), "wb") as wav:
        wav.setnchannels(signal16.shape[1])
        wav.setsampwidth(2)
        wav.setframerate(sampleRate)
        # Avoid a NumPy dependency: reinterpret contiguous int16 storage as bytes.
        wav.writeframes(bytes(signal16.view(torch.uint8).flatten().tolist()))


def make_example_mixture(sampleRate: int = 8000, seconds: float = 0.5) -> torch.Tensor:
    """Create a deterministic two-source/two-channel demonstration mixture."""
    t = torch.arange(int(sampleRate * seconds), dtype=torch.float64) / sampleRate
    sources = torch.stack((torch.sin(2 * torch.pi * 440 * t), torch.sin(2 * torch.pi * 733 * t)), dim=1)
    mixing = torch.tensor([[1.0, 0.45], [0.35, 1.0]], dtype=torch.float64)
    return sources @ mixing.mT


def main() -> None:
    """Load/create a mixture, run selected FDICA, and optionally save output."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, help="optional 16-bit PCM multichannel WAV")
    parser.add_argument("--output", type=Path, help="optional separated 16-bit PCM WAV")
    parser.add_argument("--model", choices=("LAP", "TVG"), default="LAP")
    parser.add_argument("--device", default="cpu", help='for example "cpu" or "cuda"')
    args = parser.parse_args()
    if args.input:
        obsSig, sampleRate = read_pcm16_wav(args.input, args.device)
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
    )
    print(f"input={tuple(obsSig.shape)}, separated={tuple(estSig.shape)}, model={args.model}")
    if args.output:
        write_pcm16_wav(args.output, estSig, sampleRate)


if __name__ == "__main__":
    main()
