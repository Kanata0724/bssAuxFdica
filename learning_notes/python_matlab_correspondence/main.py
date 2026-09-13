from __future__ import annotations

from pathlib import Path

import torch

from .bssAuxFdica import bssAuxFdica
from .example import read_pcm16_wav, write_pcm16_wav


def main() -> None:
    print("[main] 実行開始", flush=True)
    repo = Path(__file__).resolve().parent.parent
    # [main.m:28] [dirPath, fileName] = getInputFileNames(dataNo);
    dataDir = repo / "dataset" / "dev1_female4_src_12_E2A_conv"
    names = [
        "dev1_female4_src_1_E2A_pos050130_mic2123_conv.wav",
        "dev1_female4_src_2_E2A_pos050130_mic2123_conv.wav",
    ]
    loaded = []
    # [main.m:31] for iSrc = 1:nSrc
    for index, name in enumerate(names, start=1):
        # [main.m:32] filePath = dirPath + fileName(iSrc);
        path = dataDir / name
        print(f"[main] 音源読み込み {index}/{len(names)}: {path}", flush=True)
        # [main.m:33] [srcSig(:,:,iSrc), fs] = audioread(filePath); % srcSig: sample x mic x source
        signal, sampleRate = read_pcm16_wav(path)
        print(
            f"[main] 音源読み込み完了 {index}/{len(names)}: "
            f"shape={tuple(signal.shape)}, sample_rate={sampleRate}, "
            f"peak={signal.abs().max().item():.6f}",
            flush=True,
        )
        loaded.append((signal, sampleRate))
    images = [item[0] for item in loaded]
    sampleRate = loaded[0][1]
    if any(item[1] != sampleRate for item in loaded):
        raise ValueError("all source image files must have the same sample rate")
    print("[main] 混合信号の作成開始", flush=True)
    # [main.m:37] obsSig = sum(srcSig, 3); % obsSig: sample x mic
    obsSig = images[0] + images[1]
    # [main.m:40] peakVal = max(abs(obsSig), [], "all");
    peakVal = torch.max(torch.abs(obsSig))
    # [main.m:41] if  peakVal > 1 % clipped
    if peakVal > 1:
        # [main.m:42] obsSig = 0.99 * obsSig / peakValue; % maximum value is set to 0.99
        obsSig = 0.99 * obsSig / peakVal
        # [main.m:43] refSig = 0.99 * squeeze(srcSig(:, refMic, :)) / peakValue; % refSig: sample x source
        images = [0.99 * image / peakVal for image in images]
        print("[main] クリッピング防止の振幅正規化を適用", flush=True)
    print(f"[main] 混合信号の作成完了: shape={tuple(obsSig.shape)}", flush=True)
    print("[main] FDICA開始", flush=True)
    # [main.m:51] estSig = bssAuxFdica(obsSig, nSrc, ...
    # [main.m:52] "fftSize", fftSize, "shiftSize", shiftSize, "nIter", nIter, ...
    # [main.m:53] "isWhiten", isWhiten, "srcModel", srcModel, "refMic", refMic, ...
    # [main.m:54] "permSolver", permSolver, "isDraw", isDraw, "sampFreq", fs, "isFilt", isFilt);
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
        isDraw=True,
        sampFreq=sampleRate,
        seed=1,
        verbose=True,
    )
    print("[main] FDICA完了", flush=True)
    # [main.m:85] outDir = "./output/";
    output = Path(__file__).resolve().parent / "output"
    # [main.m:86] if ~isfolder(outDir); mkdir(outDir); end
    output.mkdir(exist_ok=True)
    # [main.m:87] audiowrite(outDir+sprintf("data%d", dataNo)+"_obs.wav", obsSig, fs); % observed signal
    # [main.m:88] audiowrite(outDir+sprintf("data%d", dataNo)+"_src1.wav", refSig(:, 1), fs); % source signal
    # [main.m:89] audiowrite(outDir+sprintf("data%d", dataNo)+"_src2.wav", refSig(:, 2), fs); % source signal
    # [main.m:90] audiowrite(outDir+sprintf("data%d", dataNo)+"_est1.wav", estSig(:, 1), fs); % estimated signal
    # [main.m:91] audiowrite(outDir+sprintf("data%d", dataNo)+"_est2.wav", estSig(:, 2), fs); % estimated signal
    outputs = [
        ("data1_obs.wav", obsSig),
        ("data1_src1.wav", images[0][:, :1]),
        ("data1_src2.wav", images[1][:, :1]),
        ("data1_est1.wav", estSig[:, :1]),
        ("data1_est2.wav", estSig[:, 1:2]),
    ]
    print(f"[main] WAV保存開始: output={output}", flush=True)
    for index, (name, signal) in enumerate(outputs, start=1):
        # [main.m:87] audiowrite(outDir+sprintf("data%d", dataNo)+"_obs.wav", obsSig, fs); % observed signal
        write_pcm16_wav(output / name, signal, sampleRate)
        print(
            f"[main] WAV保存 {index}/{len(outputs)} 完了: {name}, "
            f"peak={signal.abs().max().item():.6f}",
            flush=True,
        )
    print(f"[main] 実行完了: {len(outputs)}ファイルを保存", flush=True)


if __name__ == "__main__":
    main()
