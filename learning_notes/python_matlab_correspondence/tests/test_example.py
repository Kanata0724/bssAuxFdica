import torch
import torchaudio

from python_fdica.example import read_audio, read_pcm16_wav, write_pcm16_wav


def test_pcm16_wav_round_trip(tmp_path) -> None:
    path = tmp_path / "audio.wav"
    signal = torch.tensor([[-1.0, 0.25], [0.5, -0.5], [1.0, 0.0]], dtype=torch.float64)
    write_pcm16_wav(path, signal, 8000)
    restored, sample_rate = read_pcm16_wav(path)
    assert restored.shape == signal.shape
    assert sample_rate == 8000
    assert torch.isfinite(restored).all()


def test_read_audio_decodes_as_float32_and_restores_default_dtype(monkeypatch) -> None:
    previous = torch.get_default_dtype()
    torch.set_default_dtype(torch.float64)

    def fake_load(path: str) -> tuple[torch.Tensor, int]:
        assert path == "audio.wav"
        assert torch.get_default_dtype() == torch.float32
        return torch.tensor([[0.25, -0.5]], dtype=torch.float32), 16000

    monkeypatch.setattr(torchaudio, "load", fake_load)
    try:
        restored, sample_rate = read_audio("audio.wav")
        assert torch.get_default_dtype() == torch.float64
    finally:
        torch.set_default_dtype(previous)

    assert restored.dtype == torch.float64
    assert restored.shape == (2, 1)
    assert sample_rate == 16000
