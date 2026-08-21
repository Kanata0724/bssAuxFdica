% ======================================================================
% 学習用・日本語解説付きコピー
% 元ファイル: permSolverDoa.m
% 注意: 元コードは変更していません。空行と元コメントには追加解説を付けていません。
% 概要: 周波数ごとの分離順序を、到来方向の推定結果に基づいて揃えます。
% ======================================================================

% [解説 L1] この関数は入力を受け取り、est, perm を計算して返します。
% [意図] ファイル全体では「周波数ごとの分離順序を、到来方向の推定結果に基づいて揃えます。」という役割を担当します。
function [est, perm] = permSolverDoa(demixMat, mix, micPos, sampFreq)
% permSolverDoa solves frequency-wise permutation problem using direction
% of arrivals (DOA)
%
% [Syntax]
%  [est, perm] = permSolverDoa(mix, micPos, sampFreq)
%
% [Input]
%    demixMat: demixing matrix estimated by former FDICA (nSrc x nCh, nFreq)
%         mix: complex-valued spectrograms with permutation problem (nFreq x nTime x nSrc)
%      micPos: position of each microphone [m] (1 x nSrc)
%    sampFreq: sampling frequency of observed signal [Hz] (scalar)
%
% [Output]
%         est: permutation-aligned complex-valued estimated spectrograms (nFreq x nTime x nSrc)
%        perm: estimated permutation (nFreq x nSrc)
%
% [Note]
%    This function requires Statistics and Machine Learning Toolbox for
%    kmeans function
%
% Reference
%   H. Saruwatari, T. Kawamura, T. Nishikawa, A. Lee and K. Shikano, "Blind
%   source separation based on a fastconvergence algorithm combining ICA 
%   and beamforming," IEEE Trans. ASLP, vol. 14, no. 2, pp.666–678, 2006.
%

% Arguments check and set default values
% [解説 L29] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
arguments
% [解説 L30] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] この処理結果を同じアルゴリズム内の次段階へ渡し、計算の意味とshapeを保つためです。
    demixMat (:,:,:) {mustBeNumeric}
% [解説 L31] 処理途中で使う変数 mix を (:,:,:)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    mix (:,:,:) double
% [解説 L32] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] この処理結果を同じアルゴリズム内の次段階へ渡し、計算の意味とshapeを保つためです。
    micPos (1,:) {mustBeNonnegative}
% [解説 L33] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] この処理結果を同じアルゴリズム内の次段階へ渡し、計算の意味とshapeを保つためです。
    sampFreq (1,1) {mustBePositive}
end
% [解説 L35] この関数の複数の計算結果を受け取り、nFreq、nTime、nSrc にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[nFreq, nTime, nSrc] = size(mix, [1, 2, 3]);
% [解説 L36] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 nCh」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
nCh = size(demixMat, 2);
% [解説 L37] 条件「isreal(mix); error("'mix' must be complex-valued spectrograms.\n"); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if isreal(mix); error("'mix' must be complex-valued spectrograms.\n"); end
% [解説 L38] 条件「nCh ~= nSrc; error("nCh must be equal to nSrc.\n"); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if nCh ~= nSrc; error("nCh must be equal to nSrc.\n"); end
% [解説 L39] 条件「numel(micPos) ~= nSrc; error("numel(micPos) must be equal to size(mix, 3).\n"); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if numel(micPos) ~= nSrc; error("numel(micPos) must be equal to size(mix, 3).\n"); end
% [解説 L40] 条件「nSrc >= 3; error("permSolverDoa for nSrc >= 3 has not been implemented yet.\n"); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if nSrc >= 3; error("permSolverDoa for nSrc >= 3 has not been implemented yet.\n"); end

% Calculate DOA
% [解説 L43] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 soundSpd」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
soundSpd = 340; % sound speed [m/s]
% [解説 L44] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 freqAx」へ代入します。
% [意図] 時間信号と時間周波数表現の軸・長さ・係数を整え、DGTと逆変換を対応させるためです。
freqAx = linspace(0, sampFreq/2, nFreq).'; % nFreq x 1
% [解説 L45] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 sinDoa」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
sinDoa = zeros(nFreq, nSrc);
% [解説 L46] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 doa」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
doa = zeros(nFreq, nSrc);
% [解説 L47] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 mixMat」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
mixMat = zeros(nCh, nSrc, nFreq);
% [解説 L48] 現在処理している周波数ビン番号 を 1:nFreq の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 全周波数ビンへ同じ処理を適用し、帯域の一部だけが未処理になるのを防ぐためです。
for iFreq = 1:nFreq
% [解説 L49] 線形変換を打ち消す逆行列を求めます。結果を「処理途中で使う変数 mixMat(:, :, iFreq)」へ代入します。
% [意図] 時間信号と時間周波数表現の軸・長さ・係数を整え、DGTと逆変換を対応させるためです。
    mixMat(:, :, iFreq) = inv(demixMat(:,:,iFreq)); % mixMat(:, n, iFreq) is a steering vector of n-th source
end
% [解説 L51] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 isValid」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
isValid = zeros(nFreq, nSrc);
% [解説 L52] 現在処理している音源番号 を 1:nSrc の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
for iSrc = 1:nSrc
% [解説 L53] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 zeroCheck(1, 1)」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
    zeroCheck(1, 1) = 1;
% [解説 L54] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 zeroCheck(2:nFreq, 1)」へ代入します。
% [意図] 時間信号と時間周波数表現の軸・長さ・係数を整え、DGTと逆変換を対応させるためです。
    zeroCheck(2:nFreq, 1) = (mixMat(iSrc, 2, 2:nFreq)==0).*isinf(mixMat(iSrc, 2, 2:nFreq));
% [解説 L55] 条件を満たす要素位置を取得し、対象だけを処理します。結果を「処理途中で使う変数 nonZero」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
    nonZero = find(~zeroCheck);
% [解説 L56] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「処理途中で使う変数 srcAngle」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
    srcAngle = permute(angle(mixMat(1, iSrc, nonZero)./mixMat(2, iSrc, nonZero)), [3, 1, 2]); % nFreq x nCh x nSrc
% [解説 L57] 複素値または符号付き値を振幅・大きさへ変換します。結果を「処理途中で使う変数 sinDoa(nonZero, iSrc)」へ代入します。
% [意図] 分離品質をエネルギー比として数値化し、入力条件や手法間で公平に比較するためです。
    sinDoa(nonZero, iSrc)  = srcAngle ./ (2*pi*freqAx(nonZero) * abs(micPos(1)-micPos(2))) * soundSpd;
% [解説 L58] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 doa(nonZero, iSrc)」へ代入します。
% [意図] 分離品質をエネルギー比として数値化し、入力条件や手法間で公平に比較するためです。
    doa(nonZero, iSrc)     = real(asin(sinDoa(nonZero, iSrc))) / pi*180;
% [解説 L59] 複素値または符号付き値を振幅・大きさへ変換します。結果を「処理途中で使う変数 isValid(nonZero, iSrc)」へ代入します。
% [意図] 分離品質をエネルギー比として数値化し、入力条件や手法間で公平に比較するためです。
    isValid(nonZero, iSrc) = abs(sinDoa(nonZero, iSrc)) < 1;
end

% Apply k-means clustering
% [解説 L63] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 isValid」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
isValid = (isValid == 1); % convert to logical values
% [解説 L64] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 validDOA」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
validDOA = doa(isValid);
% [解説 L65] この関数の複数の計算結果を受け取り、centroids にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[~, centroids] = kmeans(validDOA, nSrc);
% [解説 L66] 指定した軸で平均し、ばらつきの少ない代表量を作ります。結果を「処理途中で使う変数 boundary」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
boundary = mean(centroids); % boundary DOA of clusters

% Align estimated spectrogram based on estimated permutation
% [解説 L69] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 est」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
est = zeros(nFreq, nTime, nSrc);
% [解説 L70] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「音源番号の対応を揃える置換情報」へ代入します。
% [意図] 周波数ごとに入れ替わり得る音源番号を統一し、時間波形として復元できるようにするためです。
perm = zeros(nFreq, nSrc);
% [解説 L71] 現在処理している周波数ビン番号 を 1:nFreq の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 全周波数ビンへ同じ処理を適用し、帯域の一部だけが未処理になるのを防ぐためです。
for iFreq = 1:nFreq
% [解説 L72] 条件「doa(iFreq, 1) <= boundary && doa(iFreq, 2) >= boundary」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
    if doa(iFreq, 1) <= boundary && doa(iFreq, 2) >= boundary
% [解説 L73] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「音源番号の対応を揃える置換情報」へ代入します。
% [意図] 周波数ごとに入れ替わり得る音源番号を統一し、時間波形として復元できるようにするためです。
        perm(iFreq, :) = [1, 2];
% [解説 L74] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
    else
% [解説 L75] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「音源番号の対応を揃える置換情報」へ代入します。
% [意図] 周波数ごとに入れ替わり得る音源番号を統一し、時間波形として復元できるようにするためです。
        perm(iFreq, :) = [2, 1];
    end
% [解説 L77] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 est(iFreq, :, :)」へ代入します。
% [意図] 時間信号と時間周波数表現の軸・長さ・係数を整え、DGTと逆変換を対応させるためです。
    est(iFreq, :, :) = mix(iFreq, :, perm(iFreq, :));
end

% [解説 L80] 実験の進行状況または計算結果を表示し、処理がどこまで進んだか確認できるようにします。
fprintf("Permutation solver (DOA) done.\n");
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
