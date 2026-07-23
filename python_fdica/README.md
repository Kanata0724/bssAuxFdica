# Python FDICA

このディレクトリは、リポジトリ直下のMATLAB版AuxFDICAを、関数ベースのPyTorchコードとして移植したものです。クラスは使用していません。LAP/TVGの補助関数法、周波数別白色化、projection back、COR/DOA/IPS置換解法、および任意の時間領域分離フィルタを実装しています。

## 動作環境とインストール

- Python 3.10以上（3.12で検証）
- PyTorch 2.1以上
- matplotlib 3.7以上（`isDraw=True` の描画に使用）

リポジトリのルートで次を実行します。

```bash
python -m pip install -r python_fdica/requirements.txt
```

## 入出力と軸

公開関数 `bssAuxFdica` の入力 `obsSig` は実数 `torch.Tensor` の `(sample, channel)`、出力 `estSig` は `(sample, source)` です。MATLABとの対応を優先し、`refMic` は1始まりです。内部では次の形状を使います。

| 値 | 形状 |
|---|---|
| 観測・推定スペクトログラム | `(frequency, frame, channel/source)` |
| 分離行列 | `(source, channel, frequency)` |
| oracle source image | `(sample, channel, source)` |

バッチ軸はMATLAB版に存在しないため設けていません。

## 使用方法

```python
import torch
from python_fdica import bssAuxFdica

obsSig = torch.randn(16000, 2, dtype=torch.float64)
estSig, cost = bssAuxFdica(
    obsSig,
    2,
    fftSize=1024,
    shiftSize=512,
    nIter=50,
    srcModel="LAP",       # "TVG" にするとTVG更新
    permSolver="COR",
    refMic=1,             # MATLABと同じ1始まり
    seed=1,
)
```

CPU/GPUは入力テンソルのdeviceで選択します。たとえば `obsSig = obsSig.to("cuda")` とすれば、生成される内部テンソルも同じGPUに置かれます。

最小例（生成テンソル、任意の16-bit PCM WAV入出力）は次のとおりです。

```bash
python -m python_fdica.example --model LAP --device cpu
python -m python_fdica.example --input mixture.wav --output separated.wav --model TVG
```

MATLAB `main.m` のdataset 1相当は次で実行します。

```bash
python -m python_fdica.main
```

## MATLABファイルとの対応

| MATLAB | Python | 役割 |
|---|---|---|
| `main.m` | `main.py` | dataset 1の混合、FDICA実行、保存 |
| `bssAuxFdica.m` | `bssAuxFdica.py` | 白色化、AuxFDICA、projection back、逆変換 |
| `permSolverCor.m` | `permSolverCor.py` | 大域・局所相関による置換解法 |
| `permSolverDoa.m` | `permSolverDoa.py` | 2音源のDOA置換解法 |
| `permSolverIps.m` | `permSolverIps.py` | oracle置換解法 |
| `DGTtool.m` | `stft.py` | FDICAが利用するDGT/擬似逆DGTの置換 |
| `getInputFileNames.m` | `main.py`内のdataset 1パス | 実行例で使用する入力指定 |

## LAPとTVG

`srcModel="LAP"` は各時間周波数点で `1/max(abs(Y), 10000*eps)`、`srcModel="TVG"` は `1/max(abs(Y)**2, 10000*eps)` を重みにします。それ以外のiterative-projection更新手順は共通で、MATLAB版と同じく早期停止せず `nIter` 回更新します。

## DGTtoolの置換

MATLAB版でFDICAに使われるDGTtool処理を調べ、次を `stft.py` に再現しました。

- `DGTtool("windowName", "b", ...)`: `torch.blackman_window(periodic=True)`。DGTtoolの係数定義と一致するperiodic Blackman
- `F.DGT`: hop倍数までの末尾ゼロ埋め、`windowLength-hop` サンプルの周期境界、`torch.stft(center=False, onesided=True, normalized=False)`
- `F.pinv`: フレーム列を周期拡張し、`torch.istft` の窓二乗和正規化でcanonical-dual相当の周期OLAを実行
- `isFilt=true`: 片側分離行列から共役対称な全周波数を作り、IFFT、循環シフト、線形畳み込みを実行

Python版はDGTtoolへ依存しません。

## テスト

```bash
pytest -q python_fdica/tests
```

STFT/ISTFT再構成、LAP/TVG実行、有限値、長さ、seed再現性、不正引数、IPS置換を検査します。

## MATLAB版との既知の差異・仮定

- MATLABの `local_whitening` は `nSrc < nCh` のとき固有値行列の次元が合わない記述です。Python版は意図された上位 `nSrc` 固有値・固有ベクトルだけを用いるため、過決定入力でも動作します。
- `permSolverCor.m` の局所相関一時配列は局所集合の要素数で確保した後、絶対周波数indexで書き込んでおり、条件によって範囲外になります。Python版は式の意図どおり、局所集合の要素を連続して集計します。
- 相関入力が定数列の場合、MATLAB `corr` は非有限値を返し得ます。Python版は分母を機械epsilonで下限処理し、有限なゼロ相関として扱います。
- MATLAB DOA版のk-meansはtoolbox実装・乱数状態に依存します。Python版は依存追加を避け、同じ1次元2クラスタ問題を決定的k-meansで解きます。境界付近では結果が異なる可能性があります。
- MATLABは `isDraw=false` のときcostを計算せず長さ `nIter` のゼロ配列を返し、`true` のとき初期値を含む `nIter+1` 要素へ暗黙拡張します。Python版もこの戻り値挙動を維持し、スペクトログラムとcost軌跡を描画します。
- 複数の `refMic` を指定したMATLAB出力はsource imageの4次元スペクトログラムになりますが、その後の逆変換との仕様が明確でありません。Python公開関数は時間領域出力を明確にするため、現在は参照マイク1個に限定します。
- MATLAB `main.m` のクリッピング分岐にある `peakValue` は未定義（計算値は `peakVal`）です。Pythonの実行例ではこの不明な分岐を移植していません。
