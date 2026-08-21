% ======================================================================
% 学習用・日本語解説付きコピー
% 元ファイル: bssAuxFdica.m
% 注意: 元コードは変更していません。空行と元コメントには追加解説を付けていません。
% 概要: AuxFDICA本体。STFT、分離行列の反復更新、尺度復元、置換解法、ISTFTを行います。
% ======================================================================

% [解説 L1] この関数は入力を受け取り、estSig, cost を計算して返します。
% [意図] ファイル全体では「AuxFDICA本体。STFT、分離行列の反復更新、尺度復元、置換解法、ISTFTを行います。」という役割を担当します。
function [estSig, cost] = bssAuxFdica(obsSig, nSrc, args)
% bssAuxFdica: blind source separation based on frequency-domain 
%              independent component analysis with auxiliary-function
%              technique
% [Syntax]
%   Using traditional name-value pair expression:
%   [estSig, cost] = bssAuxFdica(obsSig, nSrc, "fftSize", 1024, 
%                    "shiftSize", 512, "nIter", 50, "isWhiten", true, 
%                    "srcModel", "LAP", "refMic", 1 "permSolver", "COR", 
%                    "isDraw", false, "sampFreq", 16000)
%   Using pythonic expression (later R2021a):
%   [estSig, cost] = bssAuxFdica(obsSig, nSrc, fftSize=1024, 
%                    shiftSize=512, nIter=50, isWhiten=true, 
%                    srcModel="LAP", refMic=1, permSolver="COR", 
%                    isDraw=false, sampFreq=16000)
%
% [Input]
%      obsSig: time-domain observed signal (sample x channel)
%        nSrc: number of sources (scalar)
%     fftSize: window length in STFT (scalar, default: 1024)
%   shiftSize: window shift length in STFT (scalar, default: fftSize/2)
%       nIter: number of iterations for FDICA optimization 
%              (scalar, default: 50)
%    isWhiten: apply whitening before BSS (true/false, default: true)
%    srcModel: generative model of each source 
%              ("LAP" or "TVG", default: "LAP")
%              "LAP": isotropic complex Laplace distribution
%              "TVG": isotropic time-varying complex Gaussian distribution
%      refMic: reference microphone onto which estimated spectrogram is
%              projected by projection back technique 
%              (scalar or row vector, default: 1) 
%  permSolver: type of permutation solver 
%              ("none", "COR", "DOA", or "IPS", default: "COR")
%              "none": do not apply permutation solver after FDICA
%              "COR": correlation-based permutation solver
%              "DOA": direction-of-arrivals-based permutation solver
%              "IPS": ideal permutation solver using oracle source signals 
%                     (for checking upper bound performance of FDICA)
%      isDraw: draw spectrograms and cost function behavior 
%              (true/false, default: false)
%   sampFreq*: sampling frequency of observed signal [Hz]
%              (scalar, default: 16000, used for plotting spectrogram and 
%              DOA-based permutation solver)
% isPowRatio*: use power ratio feature for clustering (true/false, 
%              false uses raw amplitude spectrogram, default: true, used 
%              for COR-based permutation solver)
%    typeCor*: type of cost function 
%              ("Gl", "Lo", or "Gl+Lo", default: "Gl+Lo", used for
%              COR-based permutation solver)
%              "Gl": use global correlation
%              "Lo": use local correlation
%              "Gl+Lo": use both global and local correlations
%  deltaFreq*: adjacent frequencies for "Lo" correlation 
%              (scalar, default: 3, 0 means adjacent cost is not used, used
%              for COR-based permutation solver)
%  ratioFreq*: harmoinc frequencies for "Lo" correlation (scalar, 
%              if ratioFreq=3, round(iFreq/2), round(iFreq/3), 2*iFreq, 
%              and 3*iFreq and their adjacent frequencies 
%              (e.g., 2*iFreq-1 and 2*iFreq+1) are considered, default: 2,
%              0 means harmonic cost is not used, used for COR-based 
%              permutation solver)
%     micPos*: position of each microphone [m] (row vector, used for 
%              DOR-based permutation solver)
%     srcSig*: oracle source image signals used for "IPS" permutation
%              solver (sample x channel x source)
%     isFilt*: apply time-domain demixing filter to avoid circular
%              convoluation
% Arguments with * are not necessary
%
% [Output]
%      estSig: estimated signal obtained by FDICA (sample x source)
%        cost: cost function values of FDICA in each iteration (nIter x 1)
%    demixMat: estimated demixing matrix (nSrc x channel x fftSize/2+1)
%
% [Note]
%    This function requires the following functions:
%        DGTtool.m (https://github.com/KoheiYatabe/DGTtool)
%        permSolverCor.m
%        permSolverDoa.m
%        permSolverIps.m
%

% Check arguments and set default values
% [解説 L84] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
arguments
% [解説 L85] 複数音源が混ざった観測時間信号 を (:,:)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    obsSig (:,:) double
% [解説 L86] 分離対象の音源数 を (1,1)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    nSrc (1,1) double {mustBeInteger, mustBePositive}
% [解説 L87] STFTの窓長とFFT点数 を (1,1)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.fftSize (1,1) double {mustBeInteger, mustBePositive} = 1024
% [解説 L88] 隣接STFTフレーム間の移動サンプル数 を (1,1)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.shiftSize (1,1) double {mustBeInteger, mustBePositive} = 512
% [解説 L89] AuxFDICAの更新反復回数 を (1,1)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.nIter (1,1) double {mustBeInteger, mustBePositive} = 50
% [解説 L90] 事前白色化を行うかどうか を (1,1)、型 logical の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.isWhiten (1,1) logical = true
% [解説 L91] 音源振幅の確率モデル を (1,1)、型 string の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.srcModel (1,1) string {mustBeMember(args.srcModel, ["LAP", "TVG"])} = "LAP"
% [解説 L92] 尺度復元の基準にするマイク番号 を (1,:)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.refMic (1,:) double {mustBeInteger, mustBePositive} = 1
% [解説 L93] 周波数間の音源順序を揃える方法 を (1,1)、型 string の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.permSolver (1,1) string {mustBeMember(args.permSolver, ["none", "COR", "DOA", "IPS"])} = "COR"
% [解説 L94] デバッグ図を表示するかどうか を (1,1)、型 logical の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.isDraw (1,1) logical = false
% [解説 L95] 音声のサンプリング周波数 を (1,1)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.sampFreq (1,1) double {mustBePositive} = 16000
% [解説 L96] 処理途中で使う変数 isPowRatio を (1,1)、型 logical の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.isPowRatio (1,1) logical = true
% [解説 L97] 処理途中で使う変数 typeCor を (1,1)、型 string の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.typeCor (1,1) string {mustBeMember(args.typeCor, ["Gl", "Lo", "Gl+Lo"])} = "Gl+Lo"
% [解説 L98] 処理途中で使う変数 deltaFreq を (1,1)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.deltaFreq (1,1) double {mustBeInteger, mustBeNonnegative} = 3
% [解説 L99] 処理途中で使う変数 ratioFreq を (1,1)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.ratioFreq (1,1) double {mustBeInteger, mustBeNonnegative} = 2
% [解説 L100] 各マイクの位置 を (1,:)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.micPos (1,:) double {mustBeNonnegative}
% [解説 L101] 正解音源または正解音像の時間信号 を (:,:,:)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.srcSig (:,:,:) double {mustBeNumeric}
% [解説 L102] 時間領域の分離フィルタを適用するかどうか を (1,1)、型 logical の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    args.isFilt (1,1) logical = false
end
% [解説 L104] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「STFTの窓長とFFT点数」へ代入します。
% [意図] 周波数分解能と時間窓長を決め、MATLAB版のSTFT条件を固定するためです。
fftSize = args.fftSize;
% [解説 L105] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「隣接STFTフレーム間の移動サンプル数」へ代入します。
% [意図] フレーム間の重なりを決め、STFTとISTFTで同じ時間配置を使うためです。
shiftSize = args.shiftSize;
% [解説 L106] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「AuxFDICAの更新反復回数」へ代入します。
% [意図] 分離行列を十分に改善させつつ、計算時間と収束確認の範囲を固定するためです。
nIter = args.nIter;
% [解説 L107] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「事前白色化を行うかどうか」へ代入します。
% [意図] FDICA前にチャネル相関を除く処理を使うか選び、比較条件を明確にするためです。
isWhiten = args.isWhiten;
% [解説 L108] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「音源振幅の確率モデル」へ代入します。
% [意図] 振幅から補助変数を作る式とコスト関数を選び、音源の統計的な仮定を指定するためです。
srcModel = args.srcModel;
% [解説 L109] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「尺度復元の基準にするマイク番号」へ代入します。
% [意図] FDICAでは決まらない振幅と位相を、指定した実マイク上の音像へ戻すためです。
refMic = args.refMic;
% [解説 L110] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「周波数間の音源順序を揃える方法」へ代入します。
% [意図] 周波数ごとに独立に分離された音源番号を、時間波形へ戻す前に同じ順序へ揃えるためです。
permSolver = args.permSolver;
% [解説 L111] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「デバッグ図を表示するかどうか」へ代入します。
% [意図] スペクトログラムとコスト推移を表示し、分離結果と収束を確認するか選ぶためです。
isDraw = args.isDraw;
% [解説 L112] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「時間領域の分離フィルタを適用するかどうか」へ代入します。
% [意図] 周波数領域処理の循環畳み込みを避ける時間領域フィルタを使うか選ぶためです。
isFilt = args.isFilt;

% Check argument errors
% [解説 L115] この関数の複数の計算結果を受け取り、sigLen、nCh にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[sigLen, nCh] = size(obsSig, [1, 2]);
% [解説 L116] 条件「nSrc > nCh; error("'nSrc' must be equal or less than size(obsSig, 2).\n"); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if nSrc > nCh; error("'nSrc' must be equal or less than size(obsSig, 2).\n"); end
% [解説 L117] 条件「fftSize < shiftSize; error("'shiftSize' must be equal or less than fftSize.\n"); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if fftSize < shiftSize; error("'shiftSize' must be equal or less than fftSize.\n"); end
% [解説 L118] 条件「numel(refMic) > nCh; error("numel(refMic) must be equal or less than size(obsSig, 2).\n"); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if numel(refMic) > nCh; error("numel(refMic) must be equal or less than size(obsSig, 2).\n"); end

% Caluculate STFT
% [解説 L121] 時間信号と時間周波数表現の間を変換します。結果を「現在処理している周波数ビン番号」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
F = DGTtool("windowName", "b", "windowLength", fftSize, "windowShift", shiftSize); % create DGTtool instance
% [解説 L122] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「観測信号の複素時間周波数表現」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
obsSpec = F.DGT(obsSig); % STFT
% [解説 L123] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「非負周波数ビン数」へ代入します。
% [意図] 時間信号と時間周波数表現の軸・長さ・係数を整え、DGTと逆変換を対応させるためです。
nFreq = size(obsSpec, 1);

% Apply whitening (decorrelation and normalization of observed signals)
% [解説 L126] 条件「isWhiten」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if isWhiten
% [解説 L127] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obsSpecInput」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
    obsSpecInput = local_whitening(obsSpec, nSrc);
% [解説 L128] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
else
% [解説 L129] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obsSpecInput」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
    obsSpecInput = obsSpec(:, :, 1:nSrc); % discard unnecessary channels
end

% Apply FDICA
% [解説 L133] local_auxFdica が、分離スペクトログラム estSpecFdica、分離行列 demixMat、反復ごとの cost を返します。
% [意図] ここが前処理済みSTFTをAuxFDICAの反復更新へ渡す中心の呼び出しで、後続の尺度補正と収束確認にも結果を使います。
[estSpecFdica, demixMat, cost] = local_auxFdica(obsSpecInput, nIter, srcModel, isDraw);

% Apply projection back technique
% [解説 L136] projection back により、分離スペクトログラムと分離行列を基準マイク上の振幅・位相へ合わせます。
% [意図] FDICAだけでは決まらない音源ごとの複素尺度を復元し、聴取可能な音量と実マイク上の音像へ戻すためです。
[estSpecFdicaFix, demixMatFix] = local_projectionBack(estSpecFdica, obsSpec(:,:,refMic), demixMat);

% Apply permutation solver
% [解説 L139] 条件「permSolver == "none"」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if permSolver == "none"
% [解説 L140] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 estSpec」へ代入します。
% [意図] 時間信号と時間周波数表現の軸・長さ・係数を整え、DGTと逆変換を対応させるためです。
    estSpec = estSpecFdicaFix;
% [解説 L141] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 estPerm」へ代入します。
% [意図] 音源特徴の類似度と候補順を保持し、周波数間の音源番号を一貫した順序へ揃えるためです。
    estPerm = repmat(1:nSrc, [nFreq, 1]);
% [解説 L142] この処理を実行して現在の計算状態を更新し、次の処理が必要とする状態へ進めます。
elseif permSolver == "COR"
% [解説 L143] この関数の複数の計算結果を受け取り、置換補正後の分離スペクトログラム estSpec、周波数ごとの音源置換 estPerm にそれぞれ保存します。
% [意図] 並べ替えた信号と、その並べ替え規則の両方を保持して分離行列にも同じ置換を適用するためです。
    [estSpec, estPerm] = permSolverCor(estSpecFdicaFix, args.isPowRatio, args.typeCor, args.deltaFreq, args.ratioFreq);
% [解説 L144] この処理を実行して現在の計算状態を更新し、次の処理が必要とする状態へ進めます。
elseif permSolver == "DOA"
% [解説 L145] この関数の複数の計算結果を受け取り、置換補正後の分離スペクトログラム estSpec、周波数ごとの音源置換 estPerm にそれぞれ保存します。
% [意図] 並べ替えた信号と、その並べ替え規則の両方を保持して分離行列にも同じ置換を適用するためです。
    [estSpec, estPerm] = permSolverDoa(demixMatFix, estSpecFdicaFix, args.micPos, args.sampFreq);
% [解説 L146] COR・DOA以外が選ばれたため、正解音源を使うIPS置換解法の準備へ進みます。
% [意図] IPSは実運用用ではなく、周波数置換を理想的に解けた場合の上限性能を確認するための評価条件です。
else % IPS
% [解説 L147] 長さ1の軸を除き、後続処理が期待するshapeへ整えます。結果を「処理途中で使う変数 srcSpect」へ代入します。
% [意図] 時間信号と時間周波数表現の軸・長さ・係数を整え、DGTと逆変換を対応させるためです。
    srcSpect = F.DGT(squeeze(args.srcSig(:, args.refMic, :)));
% [解説 L148] この関数の複数の計算結果を受け取り、置換補正後の分離スペクトログラム estSpec、周波数ごとの音源置換 estPerm にそれぞれ保存します。
% [意図] 並べ替えた信号と、その並べ替え規則の両方を保持して分離行列にも同じ置換を適用するためです。
    [estSpec, estPerm] = permSolverIps(estSpecFdicaFix, srcSpect);
end
% [解説 L150] 現在処理している周波数ビン番号 を 1:nFreq の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 全周波数ビンへ同じ処理を適用し、帯域の一部だけが未処理になるのを防ぐためです。
for iFreq = 1:nFreq
% [解説 L151] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 demixMatFix(:, :, iFreq)」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
    demixMatFix(:, :, iFreq) = demixMatFix(estPerm(iFreq, :), :, iFreq);
end

% Calculate estimated time-domain signal
% [解説 L155] 条件「isFilt」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if isFilt
    % Apply demixing filter in time domain to avoid circular convolution
% [解説 L157] 逆行列が不安定・非正方でも扱える疑似逆行列を求めます。結果を「処理途中で使う変数 obsSigInput」へ代入します。
% [意図] 各音源像を同じマイク上で重ね、実際に分離器へ入力する混合波形を用意するためです。
    obsSigInput = F.pinv(obsSpecInput); % observed signal
% [解説 L158] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「各周波数で観測を音源へ分ける分離行列」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
    W = cat(3, demixMatFix, flip(conj(demixMatFix(:, :, 2:end-1)), 3)); % produce beyond Nyquist components
% [解説 L159] 周波数成分を時間方向の信号へ戻します。結果を「処理途中で使う変数 demixFilt」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
    demixFilt = real(ifft(W, fftSize, 3)); % fftSize x nSrc x nMic
% [解説 L160] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 demixFilt」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
    demixFilt = circshift(demixFilt, fftSize/2+1, 3); % move peak to center by circular shifting
% [解説 L161] 現在処理している音源番号 を 1:nSrc の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
    for iSrc = 1:nSrc
% [解説 L162] 処理途中で使う変数 iCh を 1:nCh の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各マイクチャネルを個別に処理した後で、マルチチャネル情報を正しく統合するためです。
        for iCh = 1:nCh
% [解説 L163] 長さ1の軸を除き、後続処理が期待するshapeへ整えます。結果を「現在処理している周波数ビン番号」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
            f = squeeze(demixFilt(iSrc, iCh, :));
% [解説 L164] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 tmp(:, iCh)」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
            tmp(:, iCh) = conv(obsSigInput(:, iCh), f); % linear convolution
        end
% [解説 L166] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「FDICAが推定した分離時間信号」へ代入します。
% [意図] 最終的に評価・保存する分離結果を得るためです。
        estSig(:, iSrc) = sum(tmp, 2);
    end
% [解説 L168] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「FDICAが推定した分離時間信号」へ代入します。
% [意図] 最終的に評価・保存する分離結果を得るためです。
    estSig(1:fftSize/2+1,:) = []; % cut initial components caused by group delay (circular shifting)
% [解説 L169] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
else
    % Calculate inverse STFT
% [解説 L171] 逆行列が不安定・非正方でも扱える疑似逆行列を求めます。結果を「FDICAが推定した分離時間信号」へ代入します。
% [意図] 最終的に評価・保存する分離結果を得るためです。
    estSig = F.pinv(estSpec);
end
% [解説 L173] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「FDICAが推定した分離時間信号」へ代入します。
% [意図] 最終的に評価・保存する分離結果を得るためです。
estSig = estSig(1:sigLen, :);

% Plot spectrograms and cost function behavior
% [解説 L176] 条件「isDraw」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if isDraw
% [解説 L177] DGTtoolの描画機能で、この段階の信号またはスペクトログラムを表示します。
% [意図] 処理段階ごとの波形・スペクトログラム・コストを比較し、分離と収束の状態を確認するためです。
    F.plot(obsSig, args.sampFreq); % observed signal
% [解説 L178] DGTtoolの描画機能で、この段階の信号またはスペクトログラムを表示します。
% [意図] 処理段階ごとの波形・スペクトログラム・コストを比較し、分離と収束の状態を確認するためです。
    F.plot(F.pinv(obsSpecInput), args.sampFreq); % observed signal input to FDICA
% [解説 L179] DGTtoolの描画機能で、この段階の信号またはスペクトログラムを表示します。
% [意図] 処理段階ごとの波形・スペクトログラム・コストを比較し、分離と収束の状態を確認するためです。
    F.plot(F.pinv(estSpecFdica), args.sampFreq); % estimated spectrogram before projection-back technique
% [解説 L180] 条件「permSolver ~= "none"」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
    if permSolver ~= "none"
% [解説 L181] DGTtoolの描画機能で、この段階の信号またはスペクトログラムを表示します。
% [意図] 処理段階ごとの波形・スペクトログラム・コストを比較し、分離と収束の状態を確認するためです。
        F.plot(F.pinv(estSpecFdicaFix), args.sampFreq); % estimated spectrogram before permutation solver
    end
% [解説 L183] DGTtoolの描画機能で、この段階の信号またはスペクトログラムを表示します。
% [意図] 処理段階ごとの波形・スペクトログラム・コストを比較し、分離と収束の状態を確認するためです。
    F.plot(estSig, args.sampFreq); % estimated signal
% [解説 L184] 反復ごとの目的関数を描き、AuxFDICAの収束挙動を表示します。
% [意図] 処理段階ごとの波形・スペクトログラム・コストを比較し、分離と収束の状態を確認するためです。
    local_plotCost(cost, nIter); % cost function behavior
end
end

%% Local functions
%--------------------------------------------------------------------------
% [解説 L190] この関数は入力を受け取り、Y, dP を計算して返します。
% [意図] ファイル全体では「AuxFDICA本体。STFT、分離行列の反復更新、尺度復元、置換解法、ISTFTを行います。」という役割を担当します。
function [Y, dP] = local_whitening(X, N)
% Whitening based on frequency-wise principal component analysis
%
% [inputs]
%    X: input spectrogram (I x J x M, nFreq x nTime x nCh)
%    N: number of sources (dimensions to which X(i,:,:) is projected)
%
% [outputs]
%    Y: output matrix (I x J x N)
%

% Initialize
% [解説 L202] この関数の複数の計算結果を受け取り、I、J にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[I, J, ~] = size(X, [1, 2, 3]); % nFreq x nTime x nCh
% [解説 L203] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「分離信号の複素時間周波数表現」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
Y = zeros(I, J, N);

% Apply frequency-wise whitening
% [解説 L206] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「処理途中で使う変数 Xp」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
Xp = permute(X, [3, 2, 1]); % I x J x M -> M x J x I
% [解説 L207] 処理途中で使う変数 i を 1:I の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 全周波数ビンへ同じ処理を適用し、帯域の一部だけが未処理になるのを防ぐためです。
for i = 1:I
% [解説 L208] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 Xi」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
    Xi = Xp(:, :, i); % M x J
% [解説 L209] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「AuxFDICA更新や射影に使う共分散・相関行列」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
    V = Xi*(Xi')/J; % covariance matrix of data matrix X (K x K)
% [解説 L210] この関数の複数の計算結果を受け取り、P、D にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
    [P, D] = eig(V); % eigenvalue decomposition (V = P*D*inv(P), P includes eigenvectors and D is a diagonal matrix with eigenvalues)
% [解説 L211] この関数の複数の計算結果を受け取り、idx にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
    [~, idx] = sort(diag(D), "descend"); % sort eigenvalues in descending order
% [解説 L212] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 D」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
    D = D(idx, idx); % sorted D
% [解説 L213] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 P」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
    P = P(:, idx); % sorted P
% [解説 L214] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 dP」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
    dP = P(:, 1:N); % top-d eigenvectors
% [解説 L215] 逆行列を明示せず線形方程式を解き、数値的に安定な係数を求めます。結果を「処理途中で使う変数 Yi」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
    Yi = sqrt(D)\(dP')*Xi; % whitened vector (N x J)
% [解説 L216] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「分離信号の複素時間周波数表現」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
    Y(i, :, :) = Yi.'; % J x N
end
end

%--------------------------------------------------------------------------
% [解説 L221] この関数は入力を受け取り、Y, W, cost を計算して返します。
% [意図] ファイル全体では「AuxFDICA本体。STFT、分離行列の反復更新、尺度復元、置換解法、ISTFTを行います。」という役割を担当します。
function [Y, W, cost] = local_auxFdica(X, nIter, srcModel, isDraw)
% BSS using FDICA
%
% [inputs]
%        X: observed spectrogram (I x J x M, nFreq x nTime x nCh, nCh=nSrc)
%    nIter: number of iterations
% srcModel: generative model of each source ("LAP" or "TVG")
%   isDraw: draw cost function behavior or not
%
% [outputs]
%        Y: estimated spectrogram (I x J x N, nFreq x nTime x nSrc)
%     cost: cost function values of FDICA in each iteration (nIter x 1)
%

% Initialize
% [解説 L236] この関数の複数の計算結果を受け取り、I、J、M にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[I, J, M] = size(X, [1,2,3]); % nFreq x nTime x nCh
% [解説 L237] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 N」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
N = M; % number of sources
% [解説 L238] 分離前に音源とチャネルを同じ順序で対応させる単位行列を作ります。結果を「処理途中で使う変数 E」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
E = repmat(eye(M), [1, 1, I]);
% [解説 L239] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「各周波数で観測を音源へ分ける分離行列」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
W = E; % initial demixing matrix (N x M x I)
% [解説 L240] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「分離信号の複素時間周波数表現」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
Y = X; % initial estimated spectrogram
% [解説 L241] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「処理途中で使う変数 Xp」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
Xp = permute(X, [3, 2, 1]); % M x J x I
% [解説 L242] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 Xph」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
Xph = pagectranspose(Xp); % J x M x I, pagewise Hermitian transpose (Xp')
% [解説 L243] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「処理途中で使う変数 Yp」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
Yp = permute(Y, [3, 2, 1]); % N x J x I
% [解説 L244] 各反復の目的関数を格納する cost をゼロで事前確保します。
% [意図] isDraw=true のとき初期値と更新後の値を記録し、AuxFDICAが目的関数を下げているか確認するためです。
cost = zeros(nIter, 1);
% [解説 L245] 条件「isDraw」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if isDraw
% [解説 L246] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「最適化または分離性能を表す評価値」へ代入します。
% [意図] 反復ごとの目的関数を記録し、更新が収束しているか確認するためです。
    cost(1,1) = local_calcFdicaCost(Yp, W, srcModel, I, J);
end

% Optimize
% [解説 L250] 実験の進行状況または計算結果を表示し、処理がどこまで進んだか確認できるようにします。
fprintf("Iteration:    ");
% [解説 L251] 現在のAuxFDICA反復番号 を 1:nIter の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
for iIter = 1:nIter
% [解説 L252] 実験の進行状況または計算結果を表示し、処理がどこまで進んだか確認できるようにします。
    fprintf("\b\b\b\b%4d", iIter);
% [解説 L253] 条件「srcModel == "LAP"」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
    if srcModel == "LAP"
% [解説 L254] LAPモデルでは、各音源・フレームの周波数方向振幅 |Y| を尺度 Rp として使い、極小値だけを 10000*eps へ丸めます。
% [意図] 次に 1./Rp を重みとして使うため、無音付近でゼロ除算や過大な重みが発生するのを防ぎます。
        Rp = max(abs(Yp), 10000*eps);
% [解説 L255] TVGモデルが選ばれた場合は、LAPとは異なる二乗振幅ベースの尺度計算へ切り替えます。
    elseif srcModel == "TVG"
% [解説 L256] TVGモデルでは各時間周波数点のパワー |Y|^2 を分散 Rp とし、極小値を 10000*eps へ丸めます。
% [意図] 時変ガウス分布の分散に応じた重みを作り、無音付近でも数値計算を安定させるためです。
        Rp = max(abs(Yp).^2, 10000*eps);
    end
    
% [解説 L259] 尺度 Rp の逆数を取り、現在小さい音源成分ほど大きく効く補助関数の重み invRp を作ります。
% [意図] 非ガウス音源モデルの目的関数を、反復射影法で解ける重み付き二次形式へ置き換えるためです。
    invRp = 1./Rp; % N x J x I
% [解説 L260] 処理途中で使う変数 n を 1:N の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
    for n = 1:N
% [解説 L261] 音源 n の時間重み invRp をマイク軸へ複製し、観測STFT Xp と要素積できる D（M x J x I）を作ります。
% [意図] 各時間フレームの観測を、現在の音源 n らしさに応じて重み付けする準備です。
        D = repmat(invRp(n, :, :), [M, 1, 1]); % M x J x I
% [解説 L262] (D.*Xp)X^H/J を周波数ページごとに計算し、音源 n 用の重み付き共分散行列 Vk（M x M x I）を作ります。
% [意図] AuxFDICAの補助関数を最小化する分離ベクトル更新に、現在の音源尺度を反映した観測統計量が必要だからです。
        Vk = pagemtimes(D.*Xp, Xph)/J; % M x M x I, pagewise matrix multiplication ((D(:,:,i).*Xp(:,:,i))*Xp(:,:,i)'/J)
% [解説 L263] 各周波数で (W Vk) wn = e_n を解き、音源 n を取り出す新しい分離ベクトル wn を求めます。
% [意図] 反復射影法の更新式 wn = (W Vk)^(-1)e_n を、逆行列を明示せず線形方程式として安定に計算するためです。
        wn = pagemldivide(pagemtimes(W, Vk), E(:, n, :)); % M x 1 x I, pagewise operation ((W(:,:,i)*Vk(:,:,i)) \ E(:, n, :))
% [解説 L264] wn^H Vk wn の平方根で割り、各周波数の分離ベクトルが Vk に関して単位ノルムになるよう正規化します。
% [意図] 更新方向を保ったまま任意の尺度を除き、補助関数法の制約を満たして数値発散を防ぐためです。
        wn = wn ./ sqrt( pagemtimes(pagemtimes(wn, "ctranspose", Vk, "none"), wn) ); % M x 1 x I, pagewise operation (wn(:,:,i)/sqrt(wn(:,:,i)'*Vk(:,:,i)*wn(:,:,i)))
% [解説 L265] 列ベクトル wn を複素共役転置し、観測へ左から掛ける行ベクトル wnh（1 x M x I）へ変換します。
% [意図] 複素STFTの分離では単なる転置ではなく共役転置が必要で、Y_n = w_n^H X を正しく計算するためです。
        wnh = pagectranspose(wn); % 1 x M x I, pagewise Hermitian transpose (wn(:,:,i)')
% [解説 L266] 更新した wnh を観測 Xp へ掛け、音源 n の分離スペクトログラム Yp(n,:,:) を直ちに更新します。
% [意図] 同じ反復内で次の音源を更新するとき、すでに改善した分離結果を使う逐次更新にするためです。
        Yp(n, :, :) = pagemtimes(wnh, Xp); % 1 x J x I, pagewise matrix multiplication (wnh(:,:,i)*Xp(:,:,i))
% [解説 L267] 分離行列 W の音源 n に対応する行を、更新済みの共役転置ベクトル wnh で置き換えます。
% [意図] 次の音源更新と次反復で、最新の分離ベクトル集合を使えるようにするためです。
        W(n, :, :) = wnh;
    end
% Readable implimentation
%     for i = 1:I
%         for n = 1:N
%             rn = Rp(n, :, i); % 1 x J
%             dg = ones(M, 1)*(1./rn); % M x J
%             Vk = (dg.*Xp(:, :, i))*Xp(:, :, i)'/J; % M x M
%             wn = (W(:, :, i)*Vk) \ E(:, n, i);
%             wn = wn/sqrt((wn')*Vk*wn);
%             Yp(n, :, i) = (wn')*Xp(:, :, i);
%             W(n, :, i) = wn';
%         end
%     end
% [解説 L281] 条件「isDraw」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
    if isDraw
% [解説 L282] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「最適化または分離性能を表す評価値」へ代入します。
% [意図] 反復ごとの目的関数を記録し、更新が収束しているか確認するためです。
        cost(iIter+1, 1) = local_calcFdicaCost(Yp, W, srcModel, I, J);
    end
end
% [解説 L285] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「分離信号の複素時間周波数表現」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
Y = permute(Yp, [3, 2, 1]);
% [解説 L286] 実験の進行状況または計算結果を表示し、処理がどこまで進んだか確認できるようにします。
fprintf(" FDICA done.\n");
end

%--------------------------------------------------------------------------
% [解説 L290] この関数は入力を受け取り、costVal を計算して返します。
% [意図] ファイル全体では「AuxFDICA本体。STFT、分離行列の反復更新、尺度復元、置換解法、ISTFTを行います。」という役割を担当します。
function costVal = local_calcFdicaCost(Yp, W, srcModel, I, J)
% [解説 L291] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 detW」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
detW = zeros(I, 1);
% [解説 L292] 処理途中で使う変数 i を 1:I の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 全周波数ビンへ同じ処理を適用し、帯域の一部だけが未処理になるのを防ぐためです。
for i = 1:I
% [解説 L293] 分離行列による体積変化を表す行列式を計算します。結果を「処理途中で使う変数 detW(i,1)」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
    detW(i,1) = det(W(:,:,i));
end
% [解説 L295] 条件「srcModel == "LAP"」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if srcModel == "LAP"
% [解説 L296] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 costVal」へ代入します。
% [意図] 反復ごとの目的関数を記録し、更新が収束しているか確認するためです。
    costVal = sum(abs(Yp), "all") - 2*J*sum(log(abs(detW)));
% [解説 L297] この処理を実行して現在の計算状態を更新し、次の処理が必要とする状態へ進めます。
elseif srcModel == "TVG"
% [解説 L298] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 costVal」へ代入します。
% [意図] 反復ごとの目的関数を記録し、更新が収束しているか確認するためです。
    costVal = sum(log(max(abs(Yp).^2, eps)), "all") - 2*J*sum(log(abs(detW)));
end
end

%--------------------------------------------------------------------------
% [解説 L303] この関数は入力を受け取り、fixY, fixW を計算して返します。
% [意図] ファイル全体では「AuxFDICA本体。STFT、分離行列の反復更新、尺度復元、置換解法、ISTFTを行います。」という役割を担当します。
function [fixY, fixW] = local_projectionBack(Y, S, W)
% Projection back technique to fix frequency-wise scales of estimated
% spectrogram obtained by FDICA
%
% [inputs]
%      Y: estimated spectrograms (I x J x N, nFreq x nTime x nSrc)
%      S: reference channel of observed spectrogram (I x J x 1)
%         or observed multichannel spectrogram (I x J x M, nFreq x nTime x nMic)
%      W: estimated emixing matrix (N x N x I, nSrc x nCh x nFreq)
%
% [outputs]
%   fixY: scale-fixed estimated spectrograms (I x J x N)
%         or scale-fitted estimated source images (I x J x N x M)
%   fixW: scale-fixed demixing matrix (N x N x I)
%         or scale-fitted demixing matrix for source images (N x N x I x M)
%

% Projection back
% [解説 L321] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「処理途中で使う変数 Yp」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
Yp = permute(Y, [3, 2, 1]); % N x J x I
% [解説 L322] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「処理途中で使う変数 Sp」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
Sp = permute(S, [3, 2, 1]); % 1 x J x 1 or M x J x I
% [解説 L323] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「処理途中で使う変数 Wp」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
Wp = permute(W, [4, 1, 2, 3]); % 1 x N x N x I
% [解説 L324] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 Yph」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
Yph = pagectranspose(Yp); % J x N x I, pagewise Hermitian transpose (Yp')
% [解説 L325] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 YpYph」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
YpYph = pagemtimes(Yp, Yph); % N x N x I, pagewise matrix multiplication (Yp*Yp')
% [解説 L326] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 YphOnYpYph」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
YphOnYpYph = pagemrdivide(Yph, YpYph); % J x N x I, pagewise matrix right-division (Yp'/(Yp*Yp'))
% [解説 L327] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 A」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
A = pagemtimes(Sp, YphOnYpYph); % 1 x N x I or M x N x I, pagewise matrix multiplication (Sp * Yp'/(Yp*Yp'))
% [解説 L328] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「処理途中で使う変数 Ap」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
Ap = permute(A, [1, 2, 4, 3]); % M x N x 1 x I
% [解説 L329] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「処理途中で使う変数 Ypp」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
Ypp = permute(Yp, [4, 1, 2, 3]); % 1 x N x J x I
% [解説 L330] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 fixY」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
fixY = Ap .* Ypp; % M x N x J x I, using implicit expansion
% [解説 L331] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「処理途中で使う変数 fixY」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
fixY = permute(fixY, [4, 3, 2, 1]); % I x J x N x M
% [解説 L332] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 fixW」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
fixW = Ap .* Wp; % M x N x N x I, using implicit expansion
% [解説 L333] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「処理途中で使う変数 fixW」へ代入します。
% [意図] 前後の行列演算で同じshapeと意味を保ち、信号処理を次の段階へ進めるためです。
fixW = permute(fixW, [2, 3, 4, 1]); % N x N x I x M

% Readable implementation
% [I, J, N] = size(Y, [1, 2, 3]); % nFreq x nTime x nSrc
% M = size(S, 3); % nCh
% A = zeros(M, N, I); % frequency-wise projection matrix
% for i = 1:I
%     for m = 1:M
%         Yi = permute(Y(i, :, :), [3, 2, 1]); % I x J x N -> N x J x 1
%         A(m, :, i) = S(i, :, m)*Yi'/(Yi*Yi');
%     end
% end
% fixY = zeros(I, J, N, M); % scale-fixed estimated spectrograms
% fixW = zeros(N, N, I, M); % scale-fixed demixing matrix
% for n = 1:N
%     for m = 1:M
%         for i = 1:I
%             fixY(i, :, n, m) = A(m, n, i)*Y(i, :, n);
%             fixW(n, :, i, m) = A(m, n, i)*W(n, :, i);
%         end
%     end
% end
end

%--------------------------------------------------------------------------
% [解説 L358] この関数は引数で受け取ったデータを処理します。
% [意図] AuxFDICA本体。STFT、分離行列の反復更新、尺度復元、置換解法、ISTFTを行います。
function local_plotCost(cost, nIter)
% [解説 L359] 信号・スペクトログラム・コスト推移を図にし、分離状態や収束挙動を目で確認できるようにします。
figure; plot(0:nIter, cost);
% [解説 L360] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
set(gca, "FontSize", 12);
% [解説 L361] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
xlabel("Number of iterations"); ylabel("Value of cost function");
% [解説 L362] 図へ補助線または枠線を加え、数値位置を読みやすくします。
% [意図] この処理結果を同じアルゴリズム内の次段階へ渡し、計算の意味とshapeを保つためです。
grid on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
