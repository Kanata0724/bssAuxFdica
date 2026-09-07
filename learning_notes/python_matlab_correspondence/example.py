# 対応注: 学習用コピー：元ファイル python_fdica/example.py（learning / b934fef 時点）。
# 対応注: MATLAB main.mとは別条件のPython実行例。入出力・分離呼び出しの役割のみ対応。
"""AuxFDICAを小さな信号や任意の音声ファイルで試すための実行例。

``python -m python_fdica.example`` で実行できる。入力ファイルを省略すると
2つの正弦波から混合信号を作るため、データセットがなくても動作確認できる。
"""

from __future__ import annotations

import argparse
import wave
from pathlib import Path

import torchaudio
import torch

from .bssAuxFdica import bssAuxFdica


def read_audio(path: str | Path, device: str = "cpu") -> tuple[torch.Tensor, int]:
    """音声ファイルをFDICA用の ``(sample, channel)`` テンソルとして読む。

    torchaudioは通常 ``(channel, sample)`` で返すため、最後に軸を転置する。
    戻り値は音声テンソルとサンプリング周波数の組である。``device`` に
    ``"cuda"`` を指定すれば、読み込み後のテンソルをGPUへ移動できる。
    """
    # TorchCodec版torchaudioは、プロセス全体の既定dtypeがfloat64だと振幅を
    # 誤って極小にする場合がある。読み込み中だけfloat32へ設定し、必ず戻す。
    defaultDtype = torch.get_default_dtype()
    try:
        torch.set_default_dtype(torch.float32)
        # [main.m:33] [srcSig(:,:,iSrc), fs] = audioread(filePath); % srcSig: sample x mic x source
        # 対応注: 読み込みの役割の対応。dtypeの一時変更とdevice指定はPython独自。
        values, sampleRate = torchaudio.load(str(path))
    finally:
        torch.set_default_dtype(defaultDtype)
    # FDICA内部の数値安定性を高めるため、軸を転置してfloat64へ変換する。
    return values.transpose(0, 1).to(device=device, dtype=torch.float64), sampleRate


# [main.m:87] audiowrite(outDir+sprintf("data%d", dataNo)+"_obs.wav", obsSig, fs); % observed signal
# 対応注: 書き出しの役割の対応。以下はPythonのPCM16変換とWAVヘッダー構築であり一行単位の移植ではない。
def write_pcm16_wav(path: str | Path, signal: torch.Tensor, sampleRate: int) -> None:
    """実数音声テンソルを一般的な16-bit PCM WAVファイルとして保存する。

    入力 ``signal`` は ``(sample, channel)`` 形状で、振幅範囲はおおむね
    -1から1を想定する。勾配情報やGPU配置はファイルへ保存できないため、
    計算グラフから切り離してCPUへ移す。
    """
    # 16-bit符号付き整数の範囲へ変換する。clampはクリッピング時の桁あふれを防ぐ。
    signal16 = (signal.detach().cpu().clamp(-1, 1) * 32767).round().to(torch.int16).contiguous()
    with wave.open(str(path), "wb") as wav:
        # WAVヘッダーへチャンネル数、1サンプル2バイト、サンプリング周波数を書く。
        wav.setnchannels(signal16.shape[1])
        wav.setsampwidth(2)
        wav.setframerate(sampleRate)
        # int16のメモリをバイト列として解釈し、全サンプルをファイルへ書き込む。
        wav.writeframes(bytes(signal16.view(torch.uint8).flatten().tolist()))


# 以前の関数名を使ったコードも壊さないため、同じ関数への別名を残す。
read_pcm16_wav = read_audio


# 対応注: Python独自の動作確認用合成信号。MATLAB main.mのデータセット混合とは別条件。
def make_example_mixture(sampleRate: int = 8000, seconds: float = 0.5) -> torch.Tensor:
    """2つの正弦波を混ぜ、再現可能な2チャンネル観測信号を作る。"""
    # arangeで各サンプルの時刻[s]を作る。
    t = torch.arange(int(sampleRate * seconds), dtype=torch.float64) / sampleRate
    # 440 Hzと733 Hzを仮想的な2音源としてsource列へ並べる。
    sources = torch.stack((torch.sin(2 * torch.pi * 440 * t), torch.sin(2 * torch.pi * 733 * t)), dim=1)
    # 各音源が各マイクへ異なる割合で届く状況を、2x2混合行列で表す。
    mixing = torch.tensor([[1.0, 0.45], [0.35, 1.0]], dtype=torch.float64)
    return sources @ mixing.mT


def main() -> None:
    """コマンドライン引数を読み、混合信号の作成・分離・保存を順に行う。"""
    # argparseを使うと、--inputなどのオプションとヘルプを安全に定義できる。
    # 対応注: Python独自：コマンドライン引数による実行例。
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, help="optional audio input supported by torchaudio")
    parser.add_argument("--output", type=Path, help="optional separated 16-bit PCM WAV")
    parser.add_argument("--model", choices=("LAP", "TVG"), default="LAP")
    parser.add_argument("--device", default="cpu", help='for example "cpu" or "cuda"')
    args = parser.parse_args()
    if args.input:
        # ファイルが指定された場合は、そのマルチチャンネル音声を観測信号にする。
        obsSig, sampleRate = read_audio(args.input, args.device)
    else:
        # 未指定時はプログラム内で短いテスト用混合信号を作る。
        sampleRate = 8000
        obsSig = make_example_mixture(sampleRate).to(args.device)
    # 2音源を10回の反復で分離する。isDraw=Trueなので結果の図も表示される。
    # [main.m:51] estSig = bssAuxFdica(obsSig, nSrc, ...
    # [main.m:52] "fftSize", fftSize, "shiftSize", shiftSize, "nIter", nIter, ...
    # [main.m:53] "isWhiten", isWhiten, "srcModel", srcModel, "refMic", refMic, ...
    # [main.m:54] "permSolver", permSolver, "isDraw", isDraw, "sampFreq", fs, "isFilt", isFilt);
    # 対応注: 分離呼び出しの役割のみ対応。FFT長・反復数などはmain.mと別条件。
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
        isDraw=True,
        verbose=True,
    )
    print(f"input={tuple(obsSig.shape)}, separated={tuple(estSig.shape)}, model={args.model}")
    if args.output:
        # --outputが指定された場合だけ、分離した2音源を1つの多チャンネルWAVへ保存する。
        write_pcm16_wav(args.output, estSig, sampleRate)


if __name__ == "__main__":
    # importされたときには実行せず、モジュールとして直接起動されたときだけ動かす。
    main()
