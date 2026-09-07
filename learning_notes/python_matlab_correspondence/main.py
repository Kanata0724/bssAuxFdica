# 対応注: 学習用コピー：元ファイル python_fdica/main.py（learning / b934fef 時点）。
# 対応注: 実行経路の学習用。main.mのBSS_EVAL評価はPython版に未実装。コピーの実行は想定しない。
"""リポジトリのdatasetを使ってAuxFDICAを一通り実行するメインプログラム。

MATLAB版 ``main.m`` のdataset 1に対応する。2つのソースイメージを読み、
加算してマイク観測信号を作り、FDICAで2音源へ分離してWAVへ保存する。
ターミナルでは各処理段階と反復回数をリアルタイムに確認できる。
"""

from __future__ import annotations

from pathlib import Path

import torch

from .bssAuxFdica import bssAuxFdica
from .example import read_pcm16_wav, write_pcm16_wav


def main() -> None:
    """dataset 1の読み込み、混合、音源分離、描画、保存を順番に実行する。"""
    print("[main] 実行開始", flush=True)
    # __file__を基準にするため、どのカレントディレクトリから起動しても
    # リポジトリ直下のdatasetを正しく見つけられる。
    repo = Path(__file__).resolve().parent.parent
    # 使用するデータセットを変える場合は、このフォルダとnamesを同じ組から選ぶ。
    # [main.m:28] [dirPath, fileName] = getInputFileNames(dataNo);
    # 対応注: Pythonはdataset 1のパスを直接記述。MATLABの8組からの選択機能はここにはない。
    dataDir = repo / "dataset" / "dev1_female4_src_12_E2A_conv"
    names = [
        "dev1_female4_src_1_E2A_pos050130_mic2123_conv.wav",
        "dev1_female4_src_2_E2A_pos050130_mic2123_conv.wav",
    ]
    # 各要素を (音声テンソル, サンプリング周波数) として保存する。
    loaded = []
    # [main.m:31] for iSrc = 1:nSrc
    for index, name in enumerate(names, start=1):
        # [main.m:32] filePath = dirPath + fileName(iSrc);
        path = dataDir / name
        print(f"[main] 音源読み込み {index}/{len(names)}: {path}", flush=True)
        # [main.m:33] [srcSig(:,:,iSrc), fs] = audioread(filePath); % srcSig: sample x mic x source
        # 対応注: Pythonは各音源をリストに保存し、三次元srcSigを作らない。
        signal, sampleRate = read_pcm16_wav(path)
        print(
            f"[main] 音源読み込み完了 {index}/{len(names)}: "
            f"shape={tuple(signal.shape)}, sample_rate={sampleRate}, "
            f"peak={signal.abs().max().item():.6f}",
            flush=True,
        )
        loaded.append((signal, sampleRate))
    # 音声テンソルだけをimagesへ、共通サンプリング周波数をsampleRateへ取り出す。
    images = [item[0] for item in loaded]
    sampleRate = loaded[0][1]
    # 対応注: Python独自：全入力のサンプリング周波数が同じか検査。
    if any(item[1] != sampleRate for item in loaded):
        raise ValueError("all source image files must have the same sample rate")
    print("[main] 混合信号の作成開始", flush=True)
    # ソースイメージは各音源が2本のマイクへ届いた信号。それらを加えると、
    # 実際にマイクで観測される2チャンネル混合信号になる。
    # [main.m:37] obsSig = sum(srcSig, 3); % obsSig: sample x mic
    # 対応注: 2音源を明示的に加算。
    obsSig = images[0] + images[1]
    # [main.m:40] peakVal = max(abs(obsSig), [], "all");
    peakVal = torch.max(torch.abs(obsSig))
    # [main.m:41] if  peakVal > 1 % clipped
    if peakVal > 1:
        # 振幅が1を超える場合は、WAV保存時に音割れしないよう全信号を同じ比率で縮小する。
        # [main.m:42] obsSig = 0.99 * obsSig / peakValue; % maximum value is set to 0.99
        # 対応注: MATLAB原文はpeakValueという未定義名。PythonはpeakValを使用。
        obsSig = 0.99 * obsSig / peakVal
        # [main.m:43] refSig = 0.99 * squeeze(srcSig(:, refMic, :)) / peakValue; % refSig: sample x source
        # 対応注: Pythonは基準micだけでなく各音源像の全channelを縮小。MATLAB原文のpeakValue表記はそのまま引用。
        images = [0.99 * image / peakVal for image in images]
        print("[main] クリッピング防止の振幅正規化を適用", flush=True)
    print(f"[main] 混合信号の作成完了: shape={tuple(obsSig.shape)}", flush=True)
    print("[main] FDICA開始", flush=True)
    # isDraw=TrueでMATLAB版と同じ各段階の図を表示し、verbose=Trueで進捗を出す。
    # [main.m:51] estSig = bssAuxFdica(obsSig, nSrc, ...
    # [main.m:52] "fftSize", fftSize, "shiftSize", shiftSize, "nIter", nIter, ...
    # [main.m:53] "isWhiten", isWhiten, "srcModel", srcModel, "refMic", refMic, ...
    # [main.m:54] "permSolver", permSolver, "isDraw", isDraw, "sampFreq", fs, "isFilt", isFilt);
    # 対応注: パラメータ値はmain.m冒頭に対応。Pythonはcostも受け取って_へ捨てる。
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
    # 分離終了後、python_fdica/outputフォルダへ確認用の5ファイルを保存する。
    # [main.m:85] outDir = "./output/";
    # 対応注: Pythonの保存先はpython_fdica/output。MATLABの./outputとは違う。main.mのBSS_EVAL評価処理はPython mainにはない。
    output = Path(__file__).resolve().parent / "output"
    # [main.m:86] if ~isfolder(outDir); mkdir(outDir); end
    output.mkdir(exist_ok=True)
    # obsは混合音、srcは正解ソースイメージ、estはFDICAの推定結果である。
    # [main.m:87] audiowrite(outDir+sprintf("data%d", dataNo)+"_obs.wav", obsSig, fs); % observed signal
    # [main.m:88] audiowrite(outDir+sprintf("data%d", dataNo)+"_src1.wav", refSig(:, 1), fs); % source signal
    # [main.m:89] audiowrite(outDir+sprintf("data%d", dataNo)+"_src2.wav", refSig(:, 2), fs); % source signal
    # [main.m:90] audiowrite(outDir+sprintf("data%d", dataNo)+"_est1.wav", estSig(:, 1), fs); % estimated signal
    # [main.m:91] audiowrite(outDir+sprintf("data%d", dataNo)+"_est2.wav", estSig(:, 2), fs); % estimated signal
    # 対応注: 以下の5つのaudiowriteをリストとループにまとめる。
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
        # 対応注: 1ファイルの書き出しに対応。
        write_pcm16_wav(output / name, signal, sampleRate)
        print(
            f"[main] WAV保存 {index}/{len(outputs)} 完了: {name}, "
            f"peak={signal.abs().max().item():.6f}",
            flush=True,
        )
    print(f"[main] 実行完了: {len(outputs)}ファイルを保存", flush=True)


# 対応注: Python独自：モジュールを直接起動した場合の入口。
if __name__ == "__main__":
    # ``python -m python_fdica.main`` で起動された場合のみmainを呼び出す。
    main()
