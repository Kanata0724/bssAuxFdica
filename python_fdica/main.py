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
    dataDir = repo / "dataset" / "dev1_female4_src_12_E2A_conv"
    names = [
        "dev1_female4_src_1_E2A_pos050130_mic2123_conv.wav",
        "dev1_female4_src_2_E2A_pos050130_mic2123_conv.wav",
    ]
    # 各要素を (音声テンソル, サンプリング周波数) として保存する。
    loaded = []
    for index, name in enumerate(names, start=1):
        path = dataDir / name
        print(f"[main] 音源読み込み {index}/{len(names)}: {path}", flush=True)
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
    if any(item[1] != sampleRate for item in loaded):
        raise ValueError("all source image files must have the same sample rate")
    print("[main] 混合信号の作成開始", flush=True)
    # ソースイメージは各音源が2本のマイクへ届いた信号。それらを加えると、
    # 実際にマイクで観測される2チャンネル混合信号になる。
    obsSig = images[0] + images[1]
    peakVal = torch.max(torch.abs(obsSig))
    if peakVal > 1:
        # 振幅が1を超える場合は、WAV保存時に音割れしないよう全信号を同じ比率で縮小する。
        obsSig = 0.99 * obsSig / peakVal
        images = [0.99 * image / peakVal for image in images]
        print("[main] クリッピング防止の振幅正規化を適用", flush=True)
    print(f"[main] 混合信号の作成完了: shape={tuple(obsSig.shape)}", flush=True)
    print("[main] FDICA開始", flush=True)
    # isDraw=TrueでMATLAB版と同じ各段階の図を表示し、verbose=Trueで進捗を出す。
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
    output = Path(__file__).resolve().parent / "output"
    output.mkdir(exist_ok=True)
    # obsは混合音、srcは正解ソースイメージ、estはFDICAの推定結果である。
    outputs = [
        ("data1_obs.wav", obsSig),
        ("data1_src1.wav", images[0][:, :1]),
        ("data1_src2.wav", images[1][:, :1]),
        ("data1_est1.wav", estSig[:, :1]),
        ("data1_est2.wav", estSig[:, 1:2]),
    ]
    print(f"[main] WAV保存開始: output={output}", flush=True)
    for index, (name, signal) in enumerate(outputs, start=1):
        write_pcm16_wav(output / name, signal, sampleRate)
        print(
            f"[main] WAV保存 {index}/{len(outputs)} 完了: {name}, "
            f"peak={signal.abs().max().item():.6f}",
            flush=True,
        )
    print(f"[main] 実行完了: {len(outputs)}ファイルを保存", flush=True)


if __name__ == "__main__":
    # ``python -m python_fdica.main`` で起動された場合のみmainを呼び出す。
    main()
