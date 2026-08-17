"""音声ファイル入出力とtorchaudio向けdtype対策のテスト。"""

import torch
import torchaudio

from python_fdica.example import read_audio, read_pcm16_wav, write_pcm16_wav


def test_pcm16_wav_round_trip(tmp_path) -> None:
    """WAVへ保存して再読込しても、形状・周波数・有限性が保たれる。"""
    # pytestのtmp_pathを使い、実際のoutputフォルダを汚さず一時ファイルを作る。
    path = tmp_path / "audio.wav"
    signal = torch.tensor([[-1.0, 0.25], [0.5, -0.5], [1.0, 0.0]], dtype=torch.float64)
    # 2チャンネルの境界値を含む短い信号を8 kHzのWAVとして往復させる。
    write_pcm16_wav(path, signal, 8000)
    restored, sample_rate = read_pcm16_wav(path)
    assert restored.shape == signal.shape
    assert sample_rate == 8000
    assert torch.isfinite(restored).all()


def test_read_audio_decodes_as_float32_and_restores_default_dtype(monkeypatch) -> None:
    """既定dtypeがfloat64でも、安全に読み込み、元の設定へ戻すことを確認する。"""
    previous = torch.get_default_dtype()
    torch.set_default_dtype(torch.float64)

    def fake_load(path: str) -> tuple[torch.Tensor, int]:
        """実ファイルの代わりに使う、torchaudio.loadの小さな模擬関数。"""
        assert path == "audio.wav"
        assert torch.get_default_dtype() == torch.float32
        return torch.tensor([[0.25, -0.5]], dtype=torch.float32), 16000

    # monkeypatchにより外部コーデックへ依存せず、関数内のdtypeだけを検証する。
    monkeypatch.setattr(torchaudio, "load", fake_load)
    try:
        restored, sample_rate = read_audio("audio.wav")
        assert torch.get_default_dtype() == torch.float64
    finally:
        torch.set_default_dtype(previous)

    # FDICAへ渡す最終テンソルはfloat64かつ(sample, channel)であるべき。
    assert restored.dtype == torch.float64
    assert restored.shape == (2, 1)
    assert sample_rate == 16000
