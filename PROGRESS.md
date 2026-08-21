# Research Progress

Last updated: 2026-08-21

## 現在地

このリポジトリは、ブラインド音源分離および AuxFDICA 関連の研究・実験を継続するための作業リポジトリである。

研究PC・自宅PC・スマートフォンからの継続作業を想定し、`PROGRESS.md` は短時間で状況を把握するための引き継ぎ資料として使う。

## 現在の研究目的

- 要確認: 具体的な研究目的、比較対象、最終的に評価したい仮説はまだ明文化されていない。
- コードから確認できる範囲では、MATLAB版 AuxFDICA と Python/PyTorch 移植版を用いて、2マイク観測の音声信号から2音源を分離する実験を行う構成になっている。
- 想定される評価軸は SDR/SIR/SAR、置換問題解法、LAP/TVG 音源モデル、白色化、時間領域分離フィルタなどだが、どれを主目的にするかは要確認。

## 作業環境・Git 状況

- 通常作業ブランチ: `codex-work`
- `main` は安定版として扱い、Codexから直接通常作業やpushを行わない。
- `output/` は実験出力用で、Gitに追加しない。
- データセットは研究入力データとして扱い、明示指示なしに変更・移動・削除しない。
- commit/push は、ユーザーから明示指示がある場合のみ行う。

## これまでに確認できる作業

- GitHub repository configured
- Research PC and home PC synchronization confirmed
- Codex desktop operation confirmed
- Codex workflow test completed
- Remote access from smartphone verified
- Remote commit/push test completed
- `AGENTS.md` は日本語化され、研究PC・自宅PC・スマートフォンからの継続作業ルールが整理されている。
- `codex-work` ブランチでの作業・commit・push の動作確認が行われている。

## 現在のコード構成

### MATLAB 版

- `main.m`
  - FDICA実験のメインスクリプト。
  - `dataset/` 内の source image WAV を読み込み、2音源を混合して観測信号を作る。
  - `bssAuxFdica.m` を実行し、`bss_eval` により SDR/SIR/SAR を評価する。
  - 推定音声などを `output/` に保存する。
- `bssAuxFdica.m`
  - AuxFDICA 本体。
  - STFT、必要に応じた白色化、AuxFDICA反復更新、projection back、置換問題解法、逆変換を行う。
  - 音源モデルは `LAP` と `TVG` に対応。
  - 置換問題解法は `none`, `COR`, `DOA`, `IPS` に対応。
- `permSolverCor.m`
  - Sawada et al. 系の相関ベース置換解法。
  - power ratio または振幅特徴を用いる。
  - 大域相関 `Gl`、局所相関 `Lo`、併用 `Gl+Lo` に対応。
- `permSolverDoa.m`
  - DOAベースの置換解法。
  - 2音源向けで、MATLAB Statistics and Machine Learning Toolbox の `kmeans` に依存する。
- `permSolverIps.m`
  - oracle source signal を使う ideal permutation solver。
  - FDICAの上限性能確認用。
- `DGTtool.m`
  - STFT/DGT および疑似逆変換に使われる外部由来ツール。
- `getInputFileNames.m`
  - `dataNo` 1-8 に対応するデータセット内WAVパスを返す。
- `bss_eval/`
  - SDR/SIR/SAR などのBSS評価関数群。

### Python / PyTorch 版

- `python_fdica/README.md`
  - Python版の使い方、MATLAB版との対応、既知の差異が整理されている。
- `python_fdica/bssAuxFdica.py`
  - MATLAB版 `bssAuxFdica.m` の関数ベースPyTorch移植。
  - LAP/TVG、白色化、projection back、COR/DOA/IPS、時間領域分離フィルタに対応。
- `python_fdica/stft.py`
  - `DGTtool.m` 相当のSTFT/ISTFT処理をPyTorchで実装。
- `python_fdica/main.py`
  - MATLAB `main.m` の dataset 1 相当を実行する例。
- `python_fdica/example.py`
  - 生成テンソルまたは音声ファイルを入力にした最小実行例。
- `python_fdica/permSolverCor.py`, `permSolverDoa.py`, `permSolverIps.py`
  - MATLAB版置換解法のPython移植。
- `python_fdica/tests/`
  - STFT/ISTFT、FDICA本体、置換解法、example入出力のテストがある。
- `python_fdica/requirements.txt`
  - `torch>=2.1`, `torchaudio>=2.1`, `matplotlib>=3.7`, `pytest>=7.4`

### データ・参考資料

- `dataset/`
  - `dev1_female4` / `dev1_male4`、`E2A` / `JR2`、source 1-4 の組み合わせから、2音源・2マイクの source image WAV が配置されている。
- `references/`
  - AuxIVA/AuxICA/置換問題に関する参考論文PDFが置かれている。
- `output/`
  - 実験出力用。Git管理対象にしない。

## 重要な結果・確認事項

- 2026-08-21 に既存 `main.m` / `python_fdica/main.py` 条件で MATLAB版とPython版の再比較を実行し、具体的な数値結果を下記に記録した。
- 要確認: Python版テストを直近で実行した日時・結果は `PROGRESS.md` には記録されていない。
- コード上、Python版にはSTFT/ISTFT再構成、LAP/TVG実行、seed再現性、時間領域分離フィルタ、描画呼び出し、不正引数、IPS置換を確認するテストが存在する。
- Python版READMEには、MATLAB版との差異や既知の仮定が明記されている。

## 2026-08-21 MATLAB/Python 既存条件再比較

### 目的

前回の独自簡略条件による比較結果を取り消し、既存コードに元から記述されている通常実行条件で、MATLAB版 AuxFDICA と Python/PyTorch版の数値的一致度、分離性能差、permutation solver差を再評価した。

### 比較条件

- 条件の出典:
  - MATLAB: `main.m` のパラメータ設定、および `getInputFileNames(dataNo=1)`
  - Python: `python_fdica/main.py` の dataset 1 パス、および `bssAuxFdica` 呼び出し
- Dataset: `dataNo=1` (`dev1_female4_src_12_E2A_conv`)
- Sources / channels: 2音源、2マイク
- `fftSize=4096`, `shiftSize=2048`, `nIter=50`
- `srcModel="LAP"`, `isWhiten=false`, `refMic=1`, `seed=1`
- 通常実行の `permSolver="COR"` を中心に、比較用に `none` と `IPS` も同じ推定スペクトログラムから評価した。
- `isDraw=true` は通常実行時の描画フラグだが、比較ツールではGUI表示を避け、同じ推定条件でcostを明示的に記録した。
- 生成物: `output/fdica_existing_conditions_compare/`

### コスト関数の一致度

- MATLAB版と現行Python版:
  - 初期cost: MATLAB `70067.79056387689`, Python `70067.79056387686`
  - 最終cost: MATLAB `-2123979.3597840276`, Python `-2120393.2504892806`
  - 絶対誤差総和: `68481.73404884388`
  - 二乗絶対誤差総和: `176639576.10559973`
  - 平均絶対誤差: `1342.7790989969387`
  - RMSE: `1861.0537606715716`
  - 最大絶対誤差: `3586.1092947470024`
  - 最大相対誤差: `0.001688391781317333`
- MATLAB式に合わせた調査用Python実装（ridge正則化なし）:
  - 最終cost: `-2124004.4061640203`
  - 絶対誤差総和: `1069.6808512662828`
  - 二乗絶対誤差総和: `78071.77017359641`
  - 平均絶対誤差: `20.974134338554567`
  - RMSE: `39.125682396754364`
  - 最大絶対誤差: `187.76308601396158`
  - 最大相対誤差: `8.84100178790064e-05`
- 解釈: 初期値と1反復後付近はほぼ一致するが、既存条件の50反復では差が蓄積する。現行Python版では `local_auxFdica` のridge正則化と `pinv` fallback がMATLAB版との差として効いている可能性が高い。

### 内部変数の一致度

- STFT `obsSpec`: 最大絶対誤差 `1.59e-14`, RMSE `4.03e-16`
- `obsSpecInput`: 最大絶対誤差 `1.59e-14`, RMSE `4.03e-16`
- MATLAB式Pythonの1反復後 `Y`: 最大絶対誤差 `1.63e-12`, RMSE `1.26e-14`
- MATLAB式Pythonの最終 `Y`: 最大絶対誤差 `4.97`, RMSE `1.68e-02`
- 現行Pythonの最終 `Y`: 最大絶対誤差 `7.60`, RMSE `8.32e-02`
- projection back後:
  - 現行Python: 最大絶対誤差 `2.46e-06`, RMSE `1.99e-08`
  - MATLAB式Python: 最大絶対誤差 `1.86e-09`, RMSE `1.51e-11`
- ISTFT後 `permSolver=none` の時間波形:
  - 現行Python: 最大絶対誤差 `1.64e-09`, RMSE `3.84e-10`
  - MATLAB式Python: 最大絶対誤差 `1.24e-12`, RMSE `2.90e-13`
- 解釈: STFT段階と1反復後は丸め誤差レベルで一致。差はAuxFDICAの反復更新を重ねる過程で増える。ただしprojection backとISTFT後の時間波形では差が大きく抑えられている。

### 分離性能

既存MATLAB `bss_eval_sources` で SDR/SIR/SAR を評価した。

- 入力混合信号:
  - SDR: `[0.1247, 0.1260]`
  - SIR: `[0.1247, 0.1260]`
- MATLAB `COR`:
  - SDR: `[6.7839, 8.3600]`
  - SIR: `[9.8228, 18.7030]`
  - SAR: `[10.1961, 8.8392]`
- 現行Python `COR`:
  - SDR: `[12.2381, 11.6432]`
  - SIR: `[23.8081, 17.7496]`
  - SAR: `[12.5698, 12.9366]`
- MATLAB `IPS`:
  - SDR: `[12.7336, 12.1234]`
  - SIR: `[24.8631, 18.1478]`
  - SAR: `[13.0222, 13.4376]`
- 現行Python `IPS`:
  - SDR: `[12.7336, 12.1234]`
  - SIR: `[24.8631, 18.1478]`
  - SAR: `[13.0222, 13.4376]`
- 解釈: `IPS` ではMATLAB版とPython版の分離性能は実質一致。通常条件の `COR` ではPython版のSDRがMATLAB版より約 `+5.45 dB`, `+3.28 dB` 高く、明確な性能差がある。

### Permutation solver

- `COR`: 2049 frequency bin 中 39 bin でMATLAB/Pythonのpermutationが不一致。
- 不一致bin（0-based）: `10, 11, 20, 41, 45, 46, 71, 72, 103, 104, 105, 1735, 1744, 1809, 1810, 1858, 1859, 1861, 1863, 1882, 1883, 1894, 1905, 1910, 1914, 1916, 1929, 1934, 1951, 1956, 1963, 1968, 1973, 1984, 1993, 1994, 2012, 2032, 2034`
- MATLAB式Pythonでも同じ39 binが不一致。
- `IPS`: 現行Python、MATLAB式Pythonともに不一致bin数 `0`。
- 解釈: 通常条件での最終性能差は、AuxFDICA本体よりも `permSolverCor` の実装差に由来する可能性が高い。

### 原因候補

1. `COR` 置換解法の実装差
   - Python版READMEに、MATLAB版 `permSolverCor.m` の局所相関一時配列と絶対周波数index処理の差異が既知事項として記録されている。
   - `IPS` では性能が一致するため、通常条件での大きな `COR` 性能差の主因候補。
2. 現行Python版 `local_auxFdica` のridge正則化
   - MATLAB本体にはない `VnReg = Vn + ridge * eye` と `pinv` fallback がある。
   - costとFDICA内部値には差を生むが、`IPS` 性能は一致しているため最終性能差の主因ではなさそう。
3. 50反復による丸め誤差・更新差の蓄積
   - STFTと1反復後はほぼ一致する一方、最終 `Y` では差が増える。
   - MATLAB式Pythonでも最終costに小差が残るため、長い反復で線形代数・更新順序の差が蓄積している可能性がある。

### 次に確認すべきこと

1. `permSolverCor.m` と `permSolverCor.py` の局所相関計算を、39個の不一致binで詳細比較する。
2. MATLAB版 `permSolverCor.m` のindex処理を「論文式どおり」と「現行MATLAB実装どおり」に分けて検証する。
3. Python版に `MATLAB互換COR` と `修正版COR` を分ける必要があるか判断する。
4. 現行Python版のridge正則化を維持するか、完全互換モードを追加するか検討する。

### 追加調査: コスト差が最終波形へ与える影響

- 条件は既存比較と同じ: `dataNo=1`, `dev1_female4_src_12_E2A_conv`, 2音源/2マイク, `fftSize=4096`, `shiftSize=2048`, `nIter=50`, `srcModel=LAP`, `isWhiten=false`, `refMic=1`, `seed=1`。
- 現行Python版のcost誤差は、絶対誤差 `>1` が14反復目、`>100` が16反復目、`>1000` が25反復目から発生。MATLAB互換Python版では `>1` が19反復目、`>100` が39反復目、`>1000` は未到達。
- STFTのRMSEは `4.03e-16`、1反復後 `Y` のRMSEは現行Python `1.17e-12` / MATLAB互換Python `1.26e-14`。差はAuxFDICA反復中に増え、25反復後は現行Python `4.79e-02` / MATLAB互換Python `1.19e-02`、50反復後は現行Python `8.32e-02` / MATLAB互換Python `1.68e-02`。
- 一方で、projection back後のRMSEは現行Python `1.99e-08`、MATLAB互換Python `1.51e-11` まで低下し、`IPS` permutation後も同じ水準を維持。ISTFT後 `IPS` 波形は現行PythonでもRMSE平均 `3.84e-10`、normalized RMSE約 `1.69e-08`-`1.77e-08`、相関ほぼ `1.0`、遅延 `0` sampleで、MATLAB版と実質一致。
- `COR` はMATLAB版とPython版で39/2049 frequency binのpermutationが異なり、不一致binのprojection-back後スペクトログラムenergy比は `0.0821`。不一致bin数は約1.9%だがenergy比は約8.2%あり、無視できない。
- `COR` permutation後スペクトログラムのRMSEは現行Python/MATLAB互換Pythonとも約 `0.3076`。ISTFT後 `COR` 波形はscale alignment後でも source 1: RMSE `8.37e-03`, normalized RMSE `0.358`, 相関 `0.934`; source 2: RMSE `7.48e-03`, normalized RMSE `0.368`, 相関 `0.930`。遅延はいずれも `0` sample。
- 結論: 現行Python版のridge正則化などによるcost/internal `Y` の差は存在するが、今回条件では `IPS` の最終分離波形・SDRにはほぼ影響しない。MATLAB/Pythonの最終性能差を支配しているのは、AuxFDICA本体より `permSolverCor` の39 bin不一致である可能性が高い。
- 次に確認すべきこと: `permSolverCor.m` / `permSolverCor.py` の39不一致binについて、局所相関値、探索窓、絶対/相対frequency index、tie breakまたは最大値選択の差を小さい範囲で比較する。

## 既知の注意点

- MATLAB `main.m` のクリッピング分岐では `peakVal` を計算しているが、分岐内で `peakValue` を参照している。クリッピングが発生する入力では未定義変数になる可能性がある。
- Python版READMEによると、MATLAB版 `local_whitening`、`permSolverCor.m` の局所相関配列、定数列相関、DOA k-means、複数 `refMic`、cost戻り値などに、MATLAB版とPython版の差異または仮定がある。
- `python_fdica/__pycache__/` や `python_fdica/tests/__pycache__/` が作業ツリー内に存在する。Git管理状況は作業前に必ず `git status` で確認する。
- `PROGRESS.md` には生ログを大量に貼らず、再開に必要な要点のみ残す。

## 未解決事項・要確認

- 研究の具体目的:
  - 要確認: どの手法・条件・評価指標を中心に研究するか。
- MATLAB版とPython版の位置づけ:
  - 要確認: Python版を主実装にするのか、MATLAB版の検証用移植として扱うのか。
- ベースライン:
  - 要確認: 比較対象とする既存手法、条件、論文設定。
- 評価条件:
  - 要確認: 使用する `dataNo`、音源組、部屋条件、FFTサイズ、反復回数、白色化有無、置換解法、音源モデル。
- 再現性:
  - 要確認: 研究PC・自宅PCで同じPython/MATLAB環境と依存関係が揃っているか。
- 実験結果:
  - 要確認: 現在の代表設定でMATLAB版とPython版がどの程度一致するか。

## 次に行うこと

1. 研究目的を1-3文で明文化する。
2. まず代表設定を1つ決める。
   - Dataset: 要確認
   - Method: 要確認
   - Parameters: 要確認
   - Number of sources: コード上は主に2音源
3. Python版テストを実行し、結果を短く記録する。
   - 例: `pytest -q python_fdica/tests`
4. MATLAB版 `main.m` の代表実験を実行できる環境か確認する。
5. MATLAB版とPython版の同一条件比較を行うか決める。
6. 実験する場合は、目的・条件・コマンド・重要な数値結果・解釈・次アクションをこのファイルへ簡潔に追記する。

## 作業再開時チェックリスト

1. `AGENTS.md` を読む。
2. `PROGRESS.md` を読む。
3. 現在ブランチが `codex-work` であることを確認する。
4. `git status` で作業ツリーを確認する。
5. 既存の未commit変更がある場合は、その内容を確認してから作業する。
6. 実験出力は `output/` に保存し、Gitに追加しない。
7. commit/push はユーザーから明示指示がある場合のみ行う。

## Codex notes

- 2026-08-20: `AGENTS.md`、既存 `PROGRESS.md`、リポジトリ構成、主要MATLAB/Pythonコード、テスト構成を確認し、研究PC・自宅PC・スマートフォンから再開しやすい引き継ぎ資料として `PROGRESS.md` を整理した。
- 今回は `PROGRESS.md` 以外のファイルは変更しない。commit/push もしない。
