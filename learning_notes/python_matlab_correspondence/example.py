from __future__ import annotations

import argparse
import wave
from pathlib import Path

import torchaudio
import torch

from .bssAuxFdica import bssAuxFdica


def read_audio(path: str | Path, device: str = "cpu") -> tuple[torch.Tensor, int]:
    defaultDtype = torch.get_default_dtype()
    try:
        torch.set_default_dtype(torch.float32)
        # [main.m:33] [srcSig(:,:,iSrc), fs] = audioread(filePath); % srcSig: sample x mic x source
        values, sampleRate = torchaudio.load(str(path))
    finally:
        torch.set_default_dtype(defaultDtype)
    return values.transpose(0, 1).to(device=device, dtype=torch.float64), sampleRate


# [main.m:87] audiowrite(outDir+sprintf("data%d", dataNo)+"_obs.wav", obsSig, fs); % observed signal
def write_pcm16_wav(path: str | Path, signal: torch.Tensor, sampleRate: int) -> None:
    signal16 = (signal.detach().cpu().clamp(-1, 1) * 32767).round().to(torch.int16).contiguous()
    with wave.open(str(path), "wb") as wav:
        wav.setnchannels(signal16.shape[1])
        wav.setsampwidth(2)
        wav.setframerate(sampleRate)
        wav.writeframes(bytes(signal16.view(torch.uint8).flatten().tolist()))


read_pcm16_wav = read_audio


def make_example_mixture(sampleRate: int = 8000, seconds: float = 0.5) -> torch.Tensor:
    t = torch.arange(int(sampleRate * seconds), dtype=torch.float64) / sampleRate
    sources = torch.stack((torch.sin(2 * torch.pi * 440 * t), torch.sin(2 * torch.pi * 733 * t)), dim=1)
    mixing = torch.tensor([[1.0, 0.45], [0.35, 1.0]], dtype=torch.float64)
    return sources @ mixing.mT


def main() -> None:
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
    # [main.m:51] estSig = bssAuxFdica(obsSig, nSrc, ...
    # [main.m:52] "fftSize", fftSize, "shiftSize", shiftSize, "nIter", nIter, ...
    # [main.m:53] "isWhiten", isWhiten, "srcModel", srcModel, "refMic", refMic, ...
    # [main.m:54] "permSolver", permSolver, "isDraw", isDraw, "sampFreq", fs, "isFilt", isFilt);
    estSig, _ = bssAuxFdica(
        obsSig,
        2,
        fftSize=256,
        shiftSize=128,
        nIter=10,
        srcModel=args.model,
        permSolver="COR",
        refMic=1,
        sampFreq=sampleRate,
        seed=1,
        isDraw=True,
        verbose=True,
    )
    print(f"input={tuple(obsSig.shape)}, separated={tuple(estSig.shape)}, model={args.model}")
    if args.output:
        write_pcm16_wav(args.output, estSig, sampleRate)


if __name__ == "__main__":
    main()
