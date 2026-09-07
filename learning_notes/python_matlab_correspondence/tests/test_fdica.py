# 対応注: 学習用コピー：元ファイル python_fdica/tests/test_fdica.py（learning / b934fef 時点）。
# 対応注: このテストファイルに直接対応するMATLABテストはない。Python独自の検証コードとして保持。
"""AuxFDICA本体のモデル、再現性、描画、入力検査を確認するテスト群。"""

import importlib

import pytest
import torch

from python_fdica import bssAuxFdica


def _mixture() -> torch.Tensor:
    """全テストで共有する、短い2音源・2チャンネル混合信号を作る。"""
    # 440 Hzと710 Hzの音源を、異なる比率で2本の仮想マイクへ混ぜる。
    t = torch.arange(320, dtype=torch.float64) / 8000
    src = torch.stack((torch.sin(2 * torch.pi * 440 * t), torch.sin(2 * torch.pi * 710 * t)), dim=1)
    return src @ torch.tensor([[1.0, 0.4], [0.3, 1.0]], dtype=torch.float64).mT


@pytest.mark.parametrize("model", ["LAP", "TVG"])
def test_lap_and_tvg_run_without_nonfinite_values(model: str) -> None:
    """LAP・TVGの両音源モデルがNaNや無限大を出さず完走する。"""
    x = _mixture()
    y, cost = bssAuxFdica(
        x, 2, fftSize=32, shiftSize=16, nIter=2, isWhiten=False,
        srcModel=model, permSolver="none", seed=7,
    )
    # determined条件では、2チャンネル入力から同じ長さの2音源が返る。
    assert y.shape == x.shape
    assert cost.shape == (2,)
    assert torch.isfinite(y).all()


def test_same_seed_is_reproducible() -> None:
    """同じ入力とseedを与えると、結果が完全に再現されることを確認する。"""
    kwargs = dict(fftSize=32, shiftSize=16, nIter=2, isWhiten=False, srcModel="LAP", permSolver="none")
    first, _ = bssAuxFdica(_mixture(), 2, seed=11, **kwargs)
    second, _ = bssAuxFdica(_mixture(), 2, seed=11, **kwargs)
    torch.testing.assert_close(first, second, atol=0, rtol=0)


def test_time_domain_demixing_filter_runs() -> None:
    """isFilt=Trueの時間領域FIR分離経路も有限値の波形を返す。"""
    y, _ = bssAuxFdica(
        _mixture(), 2, fftSize=32, shiftSize=16, nIter=1,
        isWhiten=False, srcModel="LAP", permSolver="none", isFilt=True,
    )
    assert y.shape == _mixture().shape
    assert torch.isfinite(y).all()


def test_verbose_reports_stages_and_iteration_progress(capsys: pytest.CaptureFixture[str]) -> None:
    """verbose=Trueで、主要段階と各反復の進捗が標準出力へ現れる。"""
    bssAuxFdica(
        _mixture(), 2, fftSize=32, shiftSize=16, nIter=2,
        isWhiten=False, srcModel="LAP", permSolver="none", verbose=True,
    )
    # capsysはターミナル出力を捕捉し、必要なメッセージを文字列として検査できる。
    output = capsys.readouterr().out
    assert "[FDICA] 実行開始" in output
    assert "[FDICA] STFT完了" in output
    assert "[FDICA] 学習 1/2 開始" in output
    assert "[FDICA] 学習 2/2 完了" in output
    assert "[FDICA] Projection back完了" in output
    assert "[FDICA] 実行完了" in output


def test_is_draw_calculates_cost_and_requests_plots(monkeypatch: pytest.MonkeyPatch) -> None:
    """isDraw=Trueで全描画段階と初期値を含むコストが要求される。"""
    fdica_module = importlib.import_module("python_fdica.bssAuxFdica")
    spectrogramTitles = []
    costCalls = []
    showCalls = []

    def fake_plot_spectrogram(signal: torch.Tensor, sampFreq: float, fftSize: int, shiftSize: int, *, title: str) -> None:
        """実ウィンドウを開かず、描画対象と設定だけを記録する模擬関数。"""
        assert signal.ndim == 2
        assert sampFreq == 8000
        assert fftSize == 32
        assert shiftSize == 16
        spectrogramTitles.append(title)

    def fake_plot_cost(cost: torch.Tensor, nIter: int) -> None:
        """コスト描画へ渡された配列形状と反復数を記録する。"""
        costCalls.append((cost.shape, nIter))

    # matplotlibのGUIを開かずに、描画関数が呼ばれた事実を検証する。
    monkeypatch.setattr(fdica_module, "local_plotSpectrogram", fake_plot_spectrogram)
    monkeypatch.setattr(fdica_module, "local_plotCost", fake_plot_cost)
    monkeypatch.setattr(fdica_module, "local_showPlots", lambda: showCalls.append(True))

    y, cost = bssAuxFdica(
        _mixture(), 2, fftSize=32, shiftSize=16, nIter=2, isWhiten=False,
        srcModel="LAP", permSolver="none", isDraw=True, sampFreq=8000,
    )

    assert y.shape == _mixture().shape
    assert cost.shape == (3,)
    assert torch.isfinite(cost).all()
    assert spectrogramTitles == [
        "Observed signal",
        "FDICA input signal",
        "Estimated signal before projection back",
        "Estimated signal",
    ]
    assert costCalls == [((3,), 2)]
    assert showCalls == [True]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"nSrc": 3},
        {"nSrc": 2, "srcModel": "GAUSS"},
        {"nSrc": 2, "permSolver": "bad"},
        {"nSrc": 2, "nIter": 0},
        {"nSrc": 2, "refMic": 0},
    ],
)
def test_invalid_fdica_arguments(kwargs: dict) -> None:
    """代表的な不正引数がValueErrorとして拒否されることをまとめて確認する。"""
    nSrc = kwargs.pop("nSrc")
    with pytest.raises(ValueError):
        bssAuxFdica(_mixture(), nSrc, fftSize=32, shiftSize=16, **kwargs)
