# Research Progress

Last updated: 2026-08-20

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

- 要確認: 現時点の `PROGRESS.md` には、AuxFDICA実験の具体的な数値結果は記録されていない。
- 要確認: MATLAB版 `main.m` を現在環境で最後に実行した結果は記録されていない。
- 要確認: Python版テストを直近で実行した日時・結果は `PROGRESS.md` には記録されていない。
- コード上、Python版にはSTFT/ISTFT再構成、LAP/TVG実行、seed再現性、時間領域分離フィルタ、描画呼び出し、不正引数、IPS置換を確認するテストが存在する。
- Python版READMEには、MATLAB版との差異や既知の仮定が明記されている。

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
