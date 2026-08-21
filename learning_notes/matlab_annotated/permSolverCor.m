% ======================================================================
% 学習用・日本語解説付きコピー
% 元ファイル: permSolverCor.m
% 注意: 元コードは変更していません。空行と元コメントには追加解説を付けていません。
% 概要: 周波数ごとの分離順序を、振幅・パワー系列の相関に基づいて揃えます。
% ======================================================================

% [解説 L1] この関数は入力を受け取り、est, perm を計算して返します。
% [意図] ファイル全体では「周波数ごとの分離順序を、振幅・パワー系列の相関に基づいて揃えます。」という役割を担当します。
function [est, perm] = permSolverCor(mix, isPowRatio, type, deltaFreq, ratioFreq)
% permSolverCor solves frequency-wise permutation problem using power
% ratios and clustering algorithm
%
% [Syntax]
%  [est, perm] = permSolverCor(mix)
%  [est, perm] = permSolverCor(mix, isPowRatio)
%  [est, perm] = permSolverCor(mix, isPowRatio, type)
%  [est, perm] = permSolverCor(mix, isPowRatio, type, deltaFreq)
%  [est, perm] = permSolverCor(mix, isPowRatio, type, deltaFreq, ratioFreq)
%
% [Input]
%         mix: complex-valued spectrograms with permutation problem (nFreq x nTime x nSrc)
%  isPowRatio: use power ratio feature for clustering (true/false, false uses amplitude spectrogram for clustering, default: true)
%        type: type of cost function ("Gl", "Lo", or "Gl+Lo", default: "Gl+Lo")
%   deltaFreq: adjacent frequencies for "Lo" cost (default: 3, 0 means adjacent cost is not used)
%   ratioFreq: harmoinc frequencies for "Lo" cost (if ratioFreq=3, round(iFreq/2), round(iFreq/3), 2*iFreq, and 3*iFreq and their adjacent frequencies (e.g., 2*iFreq-1 and 2*iFreq+1) are considered, default: 2, 0 means harmonic cost is not used)
%
% [Output]
%         est: permutation-aligned complex-valued estimated spectrograms (nFreq x nTime x nSrc)
%        perm: estimated permutation (nFreq x nSrc)
%
% Reference
%   H. Sawada, S. Araki, and S. Makino, "Measuring dependence of bin-wise
%   separated signals for permutation alignment in frequency-domain BSS,"
%   Proc. ISCAS, pp. 3247-3250, 2007.
%

% Arguments check and set default values
% [解説 L30] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
arguments
% [解説 L31] 処理途中で使う変数 mix を (:,:,:)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    mix (:,:,:) double
% [解説 L32] 処理途中で使う変数 isPowRatio を (1,1)、型 logical の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    isPowRatio (1,1) logical = true
% [解説 L33] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] この処理結果を同じアルゴリズム内の次段階へ渡し、計算の意味とshapeを保つためです。
    type {mustBeMember(type,{'Gl', 'Lo', 'Gl+Lo'})} = "Gl+Lo"
% [解説 L34] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] この処理結果を同じアルゴリズム内の次段階へ渡し、計算の意味とshapeを保つためです。
    deltaFreq (1,1) {mustBeNonnegative} = 3
% [解説 L35] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] この処理結果を同じアルゴリズム内の次段階へ渡し、計算の意味とshapeを保つためです。
    ratioFreq (1,1) {mustBeNonnegative} = 2
end
% [解説 L37] この関数の複数の計算結果を受け取り、nFreq、nTime、nSrc にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[nFreq, nTime, nSrc] = size(mix, [1, 2, 3]);
% [解説 L38] 条件「isreal(mix); error("'mix' must be complex-valued spectrograms.\n"); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if isreal(mix); error("'mix' must be complex-valued spectrograms.\n"); end

% Calculate feature vector for clustering
% [解説 L41] 条件「isPowRatio」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if isPowRatio
% [解説 L42] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「AuxFDICA更新や射影に使う共分散・相関行列」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
    v = abs(mix).^2 ./ sum(abs(mix).^2, 3); % power ratio Eq. (14), nFreq x nTime x nSrc
% [解説 L43] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
else
% [解説 L44] 複素値または符号付き値を振幅・大きさへ変換します。結果を「AuxFDICA更新や射影に使う共分散・相関行列」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
    v = abs(mix); % amplitude spectrogram Eq. (12), nFreq x nTime x nSrc
end

% Iterative k-means clustering
% [解説 L48] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 nPerm」へ代入します。
% [意図] 音源特徴の類似度と候補順を保持し、周波数間の音源番号を一貫した順序へ揃えるためです。
nPerm = factorial(nSrc); % number of patterns in permutation
% [解説 L49] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 allPerm」へ代入します。
% [意図] 音源特徴の類似度と候補順を保持し、周波数間の音源番号を一貫した順序へ揃えるためです。
allPerm = perms((1:nSrc)); % all permutation patterns in nSrc sources case, nPerm x nSrc
% [解説 L50] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「音源番号の対応を揃える置換情報」へ代入します。
% [意図] 周波数ごとに入れ替わり得る音源番号を統一し、時間波形として復元できるようにするためです。
perm = repmat(allPerm(end,:), [nFreq, 1]); % initial permutation so that vPerm = v, nFreq x nSrc
% [解説 L51] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 vPerm」へ代入します。
% [意図] まだ置換していない特徴量を初期状態とし、反復ごとの置換候補と比較できるようにするためです。
vPerm = v; % initial permutation-fixed feature vector
% [解説 L52] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 sumRho」へ代入します。
% [意図] 音源特徴の類似度と候補順を保持し、周波数間の音源番号を一貫した順序へ揃えるためです。
sumRho = zeros(nPerm, 1); % variable for storing cost in Eq. (18)
% [解説 L53] 実験の進行状況または計算結果を表示し、処理がどこまで進んだか確認できるようにします。
fprintf("Iteration:    "); iIter = 1;
% [解説 L54] while を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 置換候補の探索を継続し、全周波数で音源番号が安定する解を得るためです。
while(true)
% [解説 L55] 実験の進行状況または計算結果を表示し、処理がどこまで進んだか確認できるようにします。
    fprintf("\b\b\b\b%4d", iIter);

    % Store current permutation
% [解説 L58] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「音源番号の対応を揃える置換情報」へ代入します。
% [意図] 周波数ごとに入れ替わり得る音源番号を統一し、時間波形として復元できるようにするためです。
    permOld = perm; % permutation of previous iteration

% [解説 L60] 条件「type == "Gl"」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
    if type == "Gl"
        % Calculate centroid vector
% [解説 L62] 指定した軸で平均し、ばらつきの少ない代表量を作ります。結果を「処理途中で使う変数 c」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
        c = squeeze(mean(vPerm, 1)); % Eq. (17), nTime x nSrc
        
% [解説 L64] 現在処理している周波数ビン番号 を 1:nFreq の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 全周波数ビンへ同じ処理を適用し、帯域の一部だけが未処理になるのを防ぐためです。
        for iFreq = 1:nFreq
            % Calculate feature vector
% [解説 L66] 長さ1の軸を除き、後続処理が期待するshapeへ整えます。結果を「処理途中で使う変数 vf」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
            vf = squeeze(v(iFreq, :, :)); % nTime x nSrc

            % Calculate correlations between iFreq and centroid
% [解説 L69] 時間変化パターンの類似度を測り、同じ音源らしさを評価します。結果を「処理途中で使う変数 rhoGl」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
            rhoGl = corr(vf, c); % correlation between feature and centroid vectors

            % Calucule cost function value
% [解説 L72] 処理途中で使う変数 iPerm を 1:nPerm % calc Eq. (18) for all permutation patterns の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
            for iPerm = 1:nPerm % calc Eq. (18) for all permutation patterns
% [解説 L73] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 sumRho(iPerm, 1)」へ代入します。
% [意図] 音源特徴の類似度と候補順を保持し、周波数間の音源番号を一貫した順序へ揃えるためです。
                sumRho(iPerm, 1) = sum(diag(rhoGl(:, allPerm(iPerm, :)))); % diagonal elements of "rho(:, allPerm(iPerm, :))" are permuted combination
            end

            % Update permutation by maximizing cost
% [解説 L77] この関数の複数の計算結果を受け取り、perm、vPerm にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
            [perm, vPerm] = local_updatePerm(sumRho, allPerm, perm, v, vPerm, iFreq);
        end

% [解説 L80] この処理を実行して現在の計算状態を更新し、次の処理が必要とする状態へ進めます。
    elseif type == "Lo"
% [解説 L81] 現在処理している周波数ビン番号 を 1:nFreq の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 全周波数ビンへ同じ処理を適用し、帯域の一部だけが未処理になるのを防ぐためです。
        for iFreq = 1:nFreq
            % Calculate feature vector
% [解説 L83] 長さ1の軸を除き、後続処理が期待するshapeへ整えます。結果を「処理途中で使う変数 vf」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
            vf = squeeze(v(iFreq, :, :)); % nTime x nSrc

            % Define set of local frequency for vg, i.e., R(f)
% [解説 L86] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 localFreqSet」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
            localFreqSet = local_produceLocalFreqSet(iFreq, nFreq, deltaFreq, ratioFreq);
            
            % Calculate correlations between iFreq and local frequency set components
% [解説 L89] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 rhoLoFreqwise」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
            rhoLoFreqwise = local_calcCorrInLocalFreqSet(localFreqSet, nSrc, vPerm, vf);

            % Calucule cost function value
% [解説 L92] 処理途中で使う変数 iPerm を 1:nPerm % calc Eq. (19) for all permutation patterns の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
            for iPerm = 1:nPerm % calc Eq. (19) for all permutation patterns
% [解説 L93] 指定した軸で平均し、ばらつきの少ない代表量を作ります。結果を「処理途中で使う変数 rhoLo」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
                rhoLo = squeeze(mean(rhoLoFreqwise(:, :, allPerm(iPerm, :)), 1)); % diagonal elements of "rho(:, :, allPerm(iPerm, :))" are permuted combination
% [解説 L94] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 sumRho(iPerm, 1)」へ代入します。
% [意図] 音源特徴の類似度と候補順を保持し、周波数間の音源番号を一貫した順序へ揃えるためです。
                sumRho(iPerm, 1) = sum(diag(rhoLo));
            end

            % Update permutation by maximizing cost
% [解説 L98] この関数の複数の計算結果を受け取り、perm、vPerm にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
            [perm, vPerm] = local_updatePerm(sumRho, allPerm, perm, v, vPerm, iFreq);
        end
        
% [解説 L101] 直前までの条件に該当しなかった場合の代替処理へ進みます。
% [意図] この処理結果を同じアルゴリズム内の次段階へ渡し、計算の意味とshapeを保つためです。
    else % type == "Gl+Lo"
        % Calculate centroid vector
% [解説 L103] 指定した軸で平均し、ばらつきの少ない代表量を作ります。結果を「処理途中で使う変数 c」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
        c = squeeze(mean(vPerm, 1)); % Eq. (17), nTime x nSrc
        
% [解説 L105] 現在処理している周波数ビン番号 を 1:nFreq の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 全周波数ビンへ同じ処理を適用し、帯域の一部だけが未処理になるのを防ぐためです。
        for iFreq = 1:nFreq
            % Calculate feature vector
% [解説 L107] 長さ1の軸を除き、後続処理が期待するshapeへ整えます。結果を「処理途中で使う変数 vf」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
            vf = squeeze(v(iFreq, :, :)); % nTime x nSrc

            % Calculate correlations between iFreq and centroid
% [解説 L110] 時間変化パターンの類似度を測り、同じ音源らしさを評価します。結果を「処理途中で使う変数 rhoGl」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
            rhoGl = corr(vf, c); % correlation between feature and centroid vectors

            % Define set of local frequency for vg, i.e., R(f)
% [解説 L113] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 localFreqSet」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
            localFreqSet = local_produceLocalFreqSet(iFreq, nFreq, deltaFreq, ratioFreq);
            
            % Calculate correlations between iFreq and local frequency set components
% [解説 L116] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 rhoLoFreqwise」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
            rhoLoFreqwise = local_calcCorrInLocalFreqSet(localFreqSet, nSrc, vPerm, vf);

            % Calucule cost function value
% [解説 L119] 処理途中で使う変数 iPerm を 1:nPerm % calc Eq. (19) for all permutation patterns の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
            for iPerm = 1:nPerm % calc Eq. (19) for all permutation patterns
% [解説 L120] 指定した軸で平均し、ばらつきの少ない代表量を作ります。結果を「処理途中で使う変数 rhoLo」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
                rhoLo = squeeze(mean(rhoLoFreqwise(:, :, allPerm(iPerm, :)), 1)); % diagonal elements of "rho(:, :, allPerm(iPerm, :))" are permuted combination
% [解説 L121] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 sumRho(iPerm, 1)」へ代入します。
% [意図] 音源特徴の類似度と候補順を保持し、周波数間の音源番号を一貫した順序へ揃えるためです。
                sumRho(iPerm, 1) = sum(diag(rhoGl(:, allPerm(iPerm, :)))) + sum(diag(rhoLo)); % Sum of global and local costs
            end

            % Update permutation by maximizing cost
% [解説 L125] この関数の複数の計算結果を受け取り、perm、vPerm にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
            [perm, vPerm] = local_updatePerm(sumRho, allPerm, perm, v, vPerm, iFreq);
        end
    end

    % Check convergence
% [解説 L130] 条件「all(permOld==perm, 'all'); break; end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
    if all(permOld==perm, 'all'); break; end
% [解説 L131] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「現在のAuxFDICA反復番号」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
    iIter = iIter + 1;
end

% Align signal based on estimated permutation
% [解説 L135] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 est」へ代入します。
% [意図] 音源特徴または候補評価を次の判定へ渡し、permutation ambiguityを解消するためです。
est = zeros(nFreq, nTime, nSrc);
% [解説 L136] 現在処理している周波数ビン番号 を 1:nFreq の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 全周波数ビンへ同じ処理を適用し、帯域の一部だけが未処理になるのを防ぐためです。
for iFreq = 1:nFreq
% [解説 L137] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 est(iFreq, :, :)」へ代入します。
% [意図] 時間信号と時間周波数表現の軸・長さ・係数を整え、DGTと逆変換を対応させるためです。
    est(iFreq, :, :) = mix(iFreq, :, perm(iFreq,:));
end

% [解説 L140] 実験の進行状況または計算結果を表示し、処理がどこまで進んだか確認できるようにします。
fprintf(" Permutation solver (COR) done.\n");
end

%% Local functions
% [解説 L144] この関数は入力を受け取り、localSet を計算して返します。
% [意図] ファイル全体では「周波数ごとの分離順序を、振幅・パワー系列の相関に基づいて揃えます。」という役割を担当します。
function localSet = local_produceLocalFreqSet(f, F, delta, ratio)
% [解説 L145] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 adjSet」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
adjSet = [f-delta:f-1, f+1:f+delta]; % set of adjacent local frequency, i.e., A(f)
% [解説 L146] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 harSet」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
harSet = [];
% [解説 L147] 処理途中で使う変数 iRatio を 2:ratio % set of harmonic local frequency, i.e., H(f) の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
for iRatio = 2:ratio % set of harmonic local frequency, i.e., H(f)
% [解説 L148] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 harSet」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
    harSet = [harSet, round(f/iRatio)-1:round(f/iRatio)+1, f*iRatio-1:f*iRatio+1];
end
% [解説 L150] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 localSet」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
localSet = unique([adjSet, harSet]); % Union of A(f) and H(f) (and sorting)
% [解説 L151] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 localSet」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
localSet = localSet(localSet>=1 & localSet<=F); % frequency index must be in the range [1:nFreq]
end

% [解説 L154] この関数は入力を受け取り、rhoFreqwise を計算して返します。
% [意図] ファイル全体では「周波数ごとの分離順序を、振幅・パワー系列の相関に基づいて揃えます。」という役割を担当します。
function rhoFreqwise = local_calcCorrInLocalFreqSet(localSet, N, vPerm, vf)
% [解説 L155] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 rhoFreqwise」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
rhoFreqwise = zeros(numel(localSet), N, N);
% [解説 L156] 現在処理している周波数ビン番号 を localSet の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
for f = localSet
% [解説 L157] 長さ1の軸を除き、後続処理が期待するshapeへ整えます。結果を「処理途中で使う変数 vg」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
    vg = squeeze(vPerm(f, :, :)); % feature vector, nTime x nSrc
% [解説 L158] 時間変化パターンの類似度を測り、同じ音源らしさを評価します。結果を「処理途中で使う変数 rhoFreqwise(f, :, :)」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
    rhoFreqwise(f, :, :) = corr(vf, vg); % correlation between feature and vg vectors
end
end

% [解説 L162] この関数は入力を受け取り、perm, vPerm を計算して返します。
% [意図] ファイル全体では「周波数ごとの分離順序を、振幅・パワー系列の相関に基づいて揃えます。」という役割を担当します。
function [perm, vPerm] = local_updatePerm(cost, allPerm, perm, v, vPerm, f)
% [解説 L163] この関数の複数の計算結果を受け取り、idx にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[~, idx] = max(cost); % find index of maximum value
% [解説 L164] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「音源番号の対応を揃える置換情報」へ代入します。
% [意図] 周波数ごとに入れ替わり得る音源番号を統一し、時間波形として復元できるようにするためです。
perm(f, :) = allPerm(idx, :); % permutation that maximizes Eq. (18)
% [解説 L165] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 vPerm(f, :, :)」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
vPerm(f, :, :) = v(f, :, perm(f, :)); % update permutation-fixed v for calculating Eq. (17)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
