% ======================================================================
% 学習用・日本語解説付きコピー
% 元ファイル: bss_eval/bss_decomp_tvfilt.m
% 注意: 元コードは変更していません。空行と元コメントには追加解説を付けていません。
% 概要: 推定信号を正解・干渉・アーティファクト成分へ分解します。
% ======================================================================

% [解説 L1] この関数は入力を受け取り、varargout を計算して返します。
% [意図] ファイル全体では「推定信号を正解・干渉・アーティファクト成分へ分解します。」という役割を担当します。
function varargout=bss_decomp_tvfilt(varargin)

% decompose an estimated source into target/interference/noise/artefacts components, assuming the admissible distortion is a time-varying filter.
%
% Usage:
%
% [s_target,e_interf[,e_noise],e_artif]=bss_decomp_tvfilt(se,index,S[,N],tvshape,tvstep,L)
%
% Input:
%   - se: row vector of length T containing the estimated source,
%   - index: points which component of S se has to be compared to,
%   - S: n x T matrix containing the original sources,
%   - N: m x T matrix containing the noise on the obseravtions (if any).
%   - tvshape : row vector of length V at most T containing the shape of the elementary 
%     allowed time variations of the filter coefficients
%   - tvstep  : hop size (in number of samples) between two consecutive
%     variations of the filter coefficients
%   - L: the number of lags
%
% Output:
%   - s_target: row vector of length T containing the target source(s)
%   contribution,
%   - e_interf: row vector of length T containing the interferences
%   contribution,
%   - e_noise: row vector of length T containing the noise contribution (if
%   any),
%   - e_artif: row vector of length T containing the artifacts
%   contribution.
%
% Developers:  - Cedric Fevotte (cf269@cam.ac.uk) - Emmanuel Vincent
% (vincent@ircam.fr) - Remi Gribonval (remi.gribonval@irisa.fr)

% [解説 L33] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「FDICAが推定した分離時間信号」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
se=varargin{1}; index=varargin{2}; S=varargin{3};
        
% [解説 L35] 値 nargin に応じて、入力形式・音源モデル・評価モードなどの処理経路を切り替えます。
switch nargin
% [解説 L36] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 5
% [解説 L37] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 N」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
        N=[]; tvshape = varargin{4}; tvstep = varargin{5}; L = varargin{6};
% [解説 L38] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 6
% [解説 L39] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 N」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
        N=varargin{4}; tvshape = varargin{5}; tvstep = varargin{6}; L = varargin{7};
% [解説 L40] 既知の選択肢に当てはまらない入力を処理し、想定外の設定を見逃さないようにします。
    otherwise
% [解説 L41] 実験の進行状況または計算結果を表示し、処理がどこまで進んだか確認できるようにします。
        disp('Wrong number of arguments.')
end

% [解説 L44] この関数の複数の計算結果を受け取り、ne、Te にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[ne,Te]=size(se);
% [解説 L45] この関数の複数の計算結果を受け取り、n、T にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[n,T]=size(S);

%%%%%%%%%% WARNINGS %%%%%%%%%%%%%
% [解説 L48] 値 isempty(N) に応じて、入力形式・音源モデル・評価モードなどの処理経路を切り替えます。
switch isempty(N)
% [解説 L49] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 1
% [解説 L50] 条件「n>T | ne>Te, disp('Watch out: signals must be in rows.'), return; end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
        if n>T | ne>Te, disp('Watch out: signals must be in rows.'), return; end        
% [解説 L51] 条件「ne~=1, disp('Watch out: se must contain only one row.'), return; end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
        if ne~=1, disp('Watch out: se must contain only one row.'), return; end
% [解説 L52] 条件「T~=Te, disp('Watch out: se and S have different lengths.'), return; end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
        if T~=Te, disp('Watch out: se and S have different lengths.'), return; end        
% [解説 L53] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 0
% [解説 L54] この関数の複数の計算結果を受け取り、m、Tm にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
        [m,Tm]=size(N);        
% [解説 L55] 条件「n>T | ne>Te | m>Tm, disp('Watch out: signals must be in rows.'), return; end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
        if n>T | ne>Te | m>Tm, disp('Watch out: signals must be in rows.'), return; end        
% [解説 L56] 条件「ne~=1, disp('Watch out: se must contain only one row.'), return; end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
        if ne~=1, disp('Watch out: se must contain only one row.'), return; end
% [解説 L57] 条件「T~=Te, disp('Watch out: S and Se have different lengths.'), return; end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
        if T~=Te, disp('Watch out: S and Se have different lengths.'), return; end        
% [解説 L58] 条件「T~=Tm, disp('Watch out: N, S and Se have different lengths.'), return; end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
        if T~=Tm, disp('Watch out: N, S and Se have different lengths.'), return; end        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the space of target source(s)
% [解説 L64] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 target_space」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
target_space = bss_make_lags(S(index,:),L); 
% Create the space of sources
% [解説 L66] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 sources_space」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
sources_space= bss_make_lags(S,L);
% Create the noise space
% [解説 L68] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 noise_space」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
noise_space  = bss_make_lags(N,L);

% [解説 L70] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 s_target」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
s_target=zeros(1,T);
% [解説 L71] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 e_interf」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
e_interf=zeros(1,T);
% [解説 L72] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 e_artif」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
e_artif=zeros(1,T);
% [解説 L73] 条件「isempty(noise_space)==0, e_noise=zeros(1,T); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if isempty(noise_space)==0, e_noise=zeros(1,T); end

%%% Target source(s) contribution %%%
% [解説 L76] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 s_target」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
s_target = bss_tvproj(se,target_space,tvshape,tvstep);

%%% Interferences contribution %%%
% [解説 L79] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 P_S_se」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
P_S_se = bss_tvproj(se,[sources_space],tvshape,tvstep);
% [解説 L80] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 e_interf」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
e_interf = P_S_se - s_target;

% [解説 L82] 値 isempty(noise_space) に応じて、入力形式・音源モデル・評価モードなどの処理経路を切り替えます。
switch isempty(noise_space)
% [解説 L83] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 1 % No noise
        %%% Artifacts contribution %%%  
% [解説 L85] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 e_artif」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
        e_artif= se - P_S_se;
        
        %%% Output %%%
% [解説 L88] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
        varargout{1}=s_target;
% [解説 L89] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
        varargout{2}=e_interf;
% [解説 L90] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
        varargout{3}=e_artif;
        
% [解説 L92] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 0 % Noise
        %%% Noise contribution %%%
% [解説 L94] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 P_SN_se」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
        P_SN_se= bss_tvproj(se,[sources_space;noise_space],tvshape,tvstep);
% [解説 L95] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 e_noise」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
        e_noise=P_SN_se-P_S_se;
        
        %%% Artifacts contribution %%%  
% [解説 L98] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 e_artif」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
        e_artif=se-P_SN_se;
        
        %%% Output %%%
% [解説 L101] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
        varargout{1}=s_target;
% [解説 L102] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
        varargout{2}=e_interf;
% [解説 L103] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
        varargout{3}=e_noise;
% [解説 L104] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
        varargout{4}=e_artif;        
end        
