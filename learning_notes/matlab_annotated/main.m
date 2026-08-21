% ======================================================================
% 学習用・日本語解説付きコピー
% 元ファイル: main.m
% 注意: 元コードは変更していません。空行と元コメントには追加解説を付けていません。
% 概要: 実験全体の入口。入力音源から観測信号を作り、AuxFDICAを実行して性能評価と保存を行います。
% ======================================================================

% Main script for frequency-domain independent component analysis (FDICA)
% Corded by D. Kitamura (d-kitamura@ieee.org) on April 23rd, 2022

% [解説 L4] 作業領域・図・コマンド表示を初期化し、前回実行の変数や図が今回の実験へ混ざらない状態にします。
clear; close all; clc;
% [解説 L5] BSS評価関数のフォルダを検索パスへ追加し、後半のSDR・SIR・SAR計算を呼び出せるようにします。
addpath("./bss_eval/");

% Set parameters
% [解説 L8] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「再現可能な乱数系列を選ぶシード」へ代入します。
% [意図] 後で乱数生成器へ渡し、同じ実験条件を再現できるようにするためです。
seed = 1; % pseudorandom seed
% [解説 L9] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「STFTの窓長とFFT点数」へ代入します。
% [意図] 周波数分解能と時間窓長を決め、MATLAB版のSTFT条件を固定するためです。
fftSize = 4096; % window length in STFT [points]
% [解説 L10] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「隣接STFTフレーム間の移動サンプル数」へ代入します。
% [意図] フレーム間の重なりを決め、STFTとISTFTで同じ時間配置を使うためです。
shiftSize = fftSize/2; % window shift length in STFT [points]
% [解説 L11] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「分離対象の音源数」へ代入します。
% [意図] 分離行列の行数と評価対象数を、今回の2音源混合に合わせるためです。
nSrc = 2; % number of sources in observed signal
% [解説 L12] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「AuxFDICAの更新反復回数」へ代入します。
% [意図] 分離行列を十分に改善させつつ、計算時間と収束確認の範囲を固定するためです。
nIter = 50; % number of iterations of FDICA
% [解説 L13] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「事前白色化を行うかどうか」へ代入します。
% [意図] FDICA前にチャネル相関を除く処理を使うか選び、比較条件を明確にするためです。
isWhiten = false; % apply whitening before FDICA or not (true/false)
% [解説 L14] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「音源振幅の確率モデル」へ代入します。
% [意図] 振幅から補助変数を作る式とコスト関数を選び、音源の統計的な仮定を指定するためです。
srcModel = "LAP"; % generative model of each source ("LAP" or "TVG")
% [解説 L15] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「尺度復元の基準にするマイク番号」へ代入します。
% [意図] FDICAでは決まらない振幅と位相を、指定した実マイク上の音像へ戻すためです。
refMic = 1; % index of reference microphone for projection back technique
% [解説 L16] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「周波数間の音源順序を揃える方法」へ代入します。
% [意図] 周波数ごとに独立に分離された音源番号を、時間波形へ戻す前に同じ順序へ揃えるためです。
permSolver = "COR"; % type of permutation solver ("none", "COR", "DOA", or "IPS")
% [解説 L17] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「デバッグ図を表示するかどうか」へ代入します。
% [意図] スペクトログラムとコスト推移を表示し、分離結果と収束を確認するか選ぶためです。
isDraw = true; % plot spectrograms and cost function behavior for debug (true/false)
% [解説 L18] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「時間領域の分離フィルタを適用するかどうか」へ代入します。
% [意図] 周波数領域処理の循環畳み込みを避ける時間領域フィルタを使うか選ぶためです。
isFilt = false; % apply time-domain demixing filter (true/false)
% [解説 L19] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「各マイクの位置」へ代入します。
% [意図] DOA置換解法がマイク間の位相差を到来方向へ変換できるよう、配列形状を与えるためです。
micPos(1) = 0; % position of the first microphone [m]
% [解説 L20] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「各マイクの位置」へ代入します。
% [意図] DOA置換解法がマイク間の位相差を到来方向へ変換できるよう、配列形状を与えるためです。
micPos(2) = 0.0566; % position of the second microphone [m]
% [解説 L21] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「使用するデータセット組を表す番号」へ代入します。
% [意図] 読み込む音源・話者・環境の組を一つの番号で再現可能に指定するためです。
dataNo = 1; % file number of input data (see getInputFileNames) (1-8)

%% Preprocessing
% Set pseudorandom seed
% [解説 L25] 乱数シードを固定し、乱数を使う処理が毎回同じ条件から始まるようにします。
rng(seed);

% Get input file names
% [解説 L28] この関数の複数の計算結果を受け取り、入力フォルダ dirPath、音源ごとのWAV名 fileName にそれぞれ保存します。
% [意図] データ番号から入力場所と二つの音源名を一度に決め、同一条件の音源像を読み込むためです。
[dirPath, fileName] = getInputFileNames(dataNo);

% Read input source image files
% [解説 L31] 現在処理している音源番号 を 1:nSrc の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
for iSrc = 1:nSrc
% [解説 L32] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「読み込むWAVの完全なパス」へ代入します。
% [意図] 共通フォルダと音源別ファイル名を結合し、現在の音源に対応するWAVを読み込むためです。
    filePath = dirPath + fileName(iSrc);
% [解説 L33] この関数の複数の計算結果を受け取り、正解音源像 srcSig（sample x mic x source）、サンプリング周波数 fs にそれぞれ保存します。
% [意図] 波形本体と標本化条件を同時に保持し、混合生成と時間周波数変換を同じ条件で行うためです。
    [srcSig(:,:,iSrc), fs] = audioread(filePath); % srcSig: sample x mic x source
end

% Mix source images
% [解説 L37] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 各音源像を同じマイク上で重ね、実際に分離器へ入力する混合波形を用意するためです。
obsSig = sum(srcSig, 3); % obsSig: sample x mic

% Check wave clipping
% [解説 L40] 候補の最大値を取り、ピークや最良候補を求めます。結果を「観測波形の最大絶対振幅」へ代入します。
% [意図] WAVの表現範囲を超えるクリッピングが起きるか、混合後の最大振幅で判定するためです。
peakVal = max(abs(obsSig), [], "all");
% [解説 L41] 条件「peakVal > 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if  peakVal > 1 % clipped
% [解説 L42] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 波形全体の比率を保ったまま最大振幅を0.99へ収め、クリッピングを避けるためです。
% [注意] 直前に計算した変数名は peakVal ですが、ここでは peakValue を参照しています。この分岐へ入る場合は要確認です。
    obsSig = 0.99 * obsSig / peakValue; % maximum value is set to 0.99
% [解説 L43] 長さ1の軸を除き、後続処理が期待するshapeへ整えます。結果を「性能評価に使う基準マイク上の正解音源」へ代入します。
% [意図] 分離結果と同じ基準マイク上の正解を作り、公平にSDR等を評価するためです。
% [注意] この行も peakValue を参照しており、peakVal との名称不一致があります。
    refSig = 0.99 * squeeze(srcSig(:, refMic, :)) / peakValue; % refSig: sample x source
% [解説 L44] 実験の進行状況または計算結果を表示し、処理がどこまで進んだか確認できるようにします。
    fprintf('Observed signal is normalized during mixture.\n');
% [解説 L45] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
else
% [解説 L46] 長さ1の軸を除き、後続処理が期待するshapeへ整えます。結果を「性能評価に使う基準マイク上の正解音源」へ代入します。
% [意図] 分離結果と同じ基準マイク上の正解を作り、公平にSDR等を評価するためです。
    refSig = squeeze(srcSig(:, refMic, :)); % refSig: sample x source
end

%% BSS based on FDICA and permutation solver
% Sample for permSolver="COR"
% [解説 L51] 観測波形へAuxFDICAを適用して音源を推定します。結果を「FDICAが推定した分離時間信号」へ代入します。
% [意図] 最終的に評価・保存する分離結果を得るためです。
estSig = bssAuxFdica(obsSig, nSrc, ...
    "fftSize", fftSize, "shiftSize", shiftSize, "nIter", nIter, ...
    "isWhiten", isWhiten, "srcModel", srcModel, "refMic", refMic, ...
    "permSolver", permSolver, "isDraw", isDraw, "sampFreq", fs, "isFilt", isFilt);

% Sample for permSolver="DOA"
% estSig = bssAuxFdica(obsSig, nSrc, ...
%     "fftSize", fftSize, "shiftSize", shiftSize, "nIter", nIter, ...
%     "isWhiten", isWhiten, "srcModel", srcModel, "refMic", refMic, ...
%     "permSolver", permSolver, "isDraw", isDraw, "sampFreq", fs, "micPos", micPos, "isFilt", isFilt);

% Sample for permSolver="IPS"
% estSig = bssAuxFdica(obsSig, nSrc, ...
%     "fftSize", fftSize, "shiftSize", shiftSize, "nIter", nIter, ...
%     "isWhiten", isWhiten, "srcModel", srcModel, "refMic", refMic, ...
%     "permSolver", permSolver, "isDraw", isDraw, "sampFreq", fs, "srcSig", srcSig, "isFilt", isFilt);

%% Evaluation of BSS performance

% Calculate input SDR and SIR
% [解説 L71] この関数の複数の計算結果を受け取り、分離前SDR inSdr、分離前SIR inSir、分離前SAR inSar にそれぞれ保存します。
% [意図] 分離前後の目的別性能を同じ評価関数で求め、改善量とアーティファクト量を比較するためです。
[inSdr, inSir, inSar] = bss_eval_sources(repmat(obsSig(:, refMic), [1, nSrc]).', refSig.');

% Calculate output SDR, SIR, and SAR
% [解説 L74] この関数の複数の計算結果を受け取り、分離後SDR outSdr、分離後SIR outSir、分離後SAR outSar にそれぞれ保存します。
% [意図] 分離前後の目的別性能を同じ評価関数で求め、改善量とアーティファクト量を比較するためです。
[outSdr, outSir, outSar] = bss_eval_sources(estSig.', refSig.');

% Display improvements of SDR and SIR and raw SAR
% [解説 L77] 現在処理している音源番号 を 1:nSrc の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
for iSrc = 1:nSrc
% [解説 L78] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 impSdr(iSrc, 1)」へ代入します。
% [意図] 分離後SDRから分離前SDRを引き、FDICAによる歪み込みの総合改善量を求めるためです。
    impSdr(iSrc, 1) = outSdr(iSrc, 1) - inSdr(iSrc, 1);
% [解説 L79] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 impSir(iSrc, 1)」へ代入します。
% [意図] 分離後SIRから分離前SIRを引き、干渉音をどれだけ抑えたかを求めるためです。
    impSir(iSrc, 1) = outSir(iSrc, 1) - inSir(iSrc, 1);
% [解説 L80] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 rawSar(iSrc, 1)」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
    rawSar(iSrc, 1) = outSar(iSrc, 1);
% [解説 L81] 実験の進行状況または計算結果を表示し、処理がどこまで進んだか確認できるようにします。
    fprintf('Source %d\n  SDRi: %.2f[dB], SIRi: %.2f[dB], SAR: %.2f[dB]\n', iSrc, impSdr(iSrc, 1), impSir(iSrc, 1), rawSar(iSrc, 1));
end

%% Output estimated wave files
% [解説 L85] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 outDir」へ代入します。
% [意図] 実験結果をoutput配下の決まった場所へ保存し、入力コードやデータセットを汚さないためです。
outDir = "./output/";
% [解説 L86] 条件「~isfolder(outDir); mkdir(outDir); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if ~isfolder(outDir); mkdir(outDir); end
% [解説 L87] この処理を実行して現在の計算状態を更新し、次の処理が必要とする状態へ進めます。
audiowrite(outDir+sprintf("data%d", dataNo)+"_obs.wav", obsSig, fs); % observed signal
% [解説 L88] この処理を実行して現在の計算状態を更新し、次の処理が必要とする状態へ進めます。
audiowrite(outDir+sprintf("data%d", dataNo)+"_src1.wav", refSig(:, 1), fs); % source signal
% [解説 L89] この処理を実行して現在の計算状態を更新し、次の処理が必要とする状態へ進めます。
audiowrite(outDir+sprintf("data%d", dataNo)+"_src2.wav", refSig(:, 2), fs); % source signal
% [解説 L90] この処理を実行して現在の計算状態を更新し、次の処理が必要とする状態へ進めます。
audiowrite(outDir+sprintf("data%d", dataNo)+"_est1.wav", estSig(:, 1), fs); % estimated signal
% [解説 L91] この処理を実行して現在の計算状態を更新し、次の処理が必要とする状態へ進めます。
audiowrite(outDir+sprintf("data%d", dataNo)+"_est2.wav", estSig(:, 2), fs); % estimated signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
