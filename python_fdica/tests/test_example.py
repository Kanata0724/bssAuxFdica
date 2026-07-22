import torch

from python_fdica.example import read_pcm16_wav, write_pcm16_wav


def test_pcm16_wav_round_trip(tmp_path) -> None:
    path = tmp_path / "audio.wav"
    signal = torch.tensor([[-1.0, 0.25], [0.5, -0.5], [1.0, 0.0]], dtype=torch.float64)
    write_pcm16_wav(path, signal, 8000)
    restored, sample_rate = read_pcm16_wav(path)
    assert restored.shape == signal.shape
    assert sample_rate == 8000
    assert torch.isfinite(restored).all()
