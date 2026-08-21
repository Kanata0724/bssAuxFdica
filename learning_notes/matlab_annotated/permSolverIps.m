% ======================================================================
% 学習用・日本語解説付きコピー
% 元ファイル: permSolverIps.m
% 注意: 元コードは変更していません。空行と元コメントには追加解説を付けていません。
% 概要: 正解音源を参照し、評価用の理想的な周波数置換を求めます。
% ======================================================================

% [解説 L1] この関数は入力を受け取り、est, perm を計算して返します。
% [意図] ファイル全体では「正解音源を参照し、評価用の理想的な周波数置換を求めます。」という役割を担当します。
function [est, perm] = permSolverIps(mix, src)
% permSolverIps solves frequency-wise permutation problem using oracle
% source signals (ideal permutation solver)
%
% [Syntax]
%  [est, perm] = permSolverIps(mix, src)
%
% [Input]
%         mix: complex-valued spectrograms with permutation problem (nFreq x nTime x nSrc)
%         src: oracle source spectrogram (nFreq x nTime x nSrc)
%
% [Output]
%         est: permutation-aligned complex-valued estimated spectrograms (nFreq x nTime x nSrc)
%        perm: estimated permutation (nFreq x nSrc)
%

% Arguments check and set default values
% [解説 L18] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
arguments
% [解説 L19] 処理途中で使う変数 mix を (:,:,:)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    mix (:,:,:) double
% [解説 L20] 処理途中で使う変数 src を (:,:,:)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    src (:,:,:) double
end
% [解説 L22] この関数の複数の計算結果を受け取り、nFreq、nTime、nSrc にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[nFreq, nTime, nSrc] = size(mix, [1, 2, 3]);
% [解説 L23] 条件「isreal(mix); error("'mix' must be complex-valued spectrograms.\n"); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if isreal(mix); error("'mix' must be complex-valued spectrograms.\n"); end
% [解説 L24] 条件「isreal(src); error("'src' must be complex-valued spectrograms.\n"); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if isreal(src); error("'src' must be complex-valued spectrograms.\n"); end
% [解説 L25] 条件「~isequal(size(mix), size(src)); error("Sizes of 'mix' and 'src' must be equal.\n"); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if ~isequal(size(mix), size(src)); error("Sizes of 'mix' and 'src' must be equal.\n"); end

% Align estimated spectrogram using oracle source spectrogram
% [解説 L28] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 est」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
est = zeros(nFreq, nTime, nSrc);
% [解説 L29] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「音源番号の対応を揃える置換情報」へ代入します。
% [意図] 周波数ごとに入れ替わり得る音源番号を統一し、時間波形として復元できるようにするためです。
perm = zeros(nFreq, nSrc);
% [解説 L30] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「音源番号の対応を揃える置換情報」へ代入します。
% [意図] 周波数ごとに入れ替わり得る音源番号を統一し、時間波形として復元できるようにするためです。
permAll = perms(1:nSrc);
% [解説 L31] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 nPerm」へ代入します。
% [意図] 音源特徴の類似度と候補順を保持し、周波数間の音源番号を一貫した順序へ揃えるためです。
nPerm = size(permAll, 1); % = factorial(nSrc)
% [解説 L32] 現在処理している周波数ビン番号 を 1:nFreq の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 全周波数ビンへ同じ処理を適用し、帯域の一部だけが未処理になるのを防ぐためです。
for iFreq = 1:nFreq
% [解説 L33] 処理途中で使う変数 iPerm を 1:nPerm の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
    for iPerm = 1:nPerm
% [解説 L34] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 err(iPerm)」へ代入します。
% [意図] 音源特徴の類似度と候補順を保持し、周波数間の音源番号を一貫した順序へ揃えるためです。
        err(iPerm) = 0;
% [解説 L35] 現在処理している音源番号 を 1:nSrc の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
        for iSrc = 1:nSrc
% [解説 L36] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 mixInd」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
            mixInd = permAll(iPerm, iSrc);
% [解説 L37] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 err(iPerm)」へ代入します。
% [意図] 音源特徴の類似度と候補順を保持し、周波数間の音源番号を一貫した順序へ揃えるためです。
            err(iPerm) = err(iPerm) + sum(abs(mix(iFreq, :, mixInd) - src(iFreq, :, iSrc)).^2, "all");
        end
    end
% [解説 L40] この関数の複数の計算結果を受け取り、ind にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
    [~, ind] = min(err);
% [解説 L41] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「音源番号の対応を揃える置換情報」へ代入します。
% [意図] 周波数ごとに入れ替わり得る音源番号を統一し、時間波形として復元できるようにするためです。
    perm(iFreq, :) = permAll(ind, :);
% [解説 L42] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 est(iFreq, :, :)」へ代入します。
% [意図] 時間信号と時間周波数表現の軸・長さ・係数を整え、DGTと逆変換を対応させるためです。
    est(iFreq, :, :) = mix(iFreq, :, perm(iFreq,:));
end

% [解説 L45] 実験の進行状況または計算結果を表示し、処理がどこまで進んだか確認できるようにします。
fprintf("Permutation solver (IPS) done.\n");
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
