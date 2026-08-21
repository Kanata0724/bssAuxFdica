% ======================================================================
% 学習用・日本語解説付きコピー
% 元ファイル: bss_eval/bss_crit.m
% 注意: 元コードは変更していません。空行と元コメントには追加解説を付けていません。
% 概要: 全候補の評価値から最適な音源対応を探索します。
% ======================================================================

% [解説 L1] この関数は入力を受け取り、varargout を計算して返します。
% [意図] ファイル全体では「全候補の評価値から最適な音源対応を探索します。」という役割を担当します。
function varargout=bss_crit(varargin)

% compute evaluation criteria given a decomposition of an estimated source into target/interference/noise/artifacts 
% of the form
%
%   se = s_target + e_interf (+ e_noise) + e_artif 
%
% Usage:
%
% 1) Global mode
%
% [SDR,SIR,(SNR,)SAR]=bss_crit(s_target,e_interf[,e_noise],e_artif)
%
% Input:
%   - s_target: row vector of length T containing the target source(s)
%   contribution,
%   - e_interf: row vector of length T containing the interferences
%   contribution,
%   - e_noise: row vector of length T containing the noise contribution 
%   (if any),
%   - e_artif: row vector of length T containing the artifacts
%   contribution.
%
% Output:
%   - SDR: Source to Distortion Ratio,
%   - SIR: Source to Interferences Ratio,
%   - SNR: Signal to Noise Ratio (if e_noise is provided),
%   - SAR: Source to Artifacts Ratio.
%
% 2) Local mode
%
% [SDR,SIR,(SNR,)SAR]=bss_crit(s_target,e_interf[,e_noise],e_artif,WINDOW,NOVERLAP)
%
% Additional input:
%   - WINDOW: 1 x W window
%   - NOVERLAP: number of samples of overlap between consecutive windows
%
% Output:
%   - SDR: n_frames x 1 vector containing local Source to Distortion Ratio,
%   - SIR: n_frames x 1 vector containing local Source to Interferences Ratio,
%   - SNR: n_frames x 1 vector containing local Signal to Noise Ratio,
%   - SAR: n_frames x 1 vector containing local Source to Artifacts Ratio.
%
% Developers:  - Cedric Fevotte (cf269@cam.ac.uk) - Emmanuel Vincent
% (vincent@ircam.fr) - Remi Gribonval (remi.gribonval@irisa.fr)


% [解説 L48] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 s_target」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
s_target=varargin{1}; e_interf=varargin{2}; 

% [解説 L50] 値 nargin に応じて、入力形式・音源モデル・評価モードなどの処理経路を切り替えます。
switch nargin
% [解説 L51] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 3
% [解説 L52] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 e_noise」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
        e_noise=[]; e_artif=varargin{3};
% [解説 L53] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 mode」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
        mode='global';
% [解説 L54] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 4
% [解説 L55] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 e_noise」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
        e_noise=varargin{3}; e_artif=varargin{4};
% [解説 L56] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 mode」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
        mode='global';
% [解説 L57] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 5
% [解説 L58] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 e_noise」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
        e_noise=[]; e_artif=varargin{3};
% [解説 L59] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 WINDOW」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
        WINDOW=varargin{4}; NOVERLAP=varargin{5};
% [解説 L60] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 mode」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
        mode='local';
% [解説 L61] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 6
% [解説 L62] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 e_noise」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
        e_noise=varargin{3}; e_artif=varargin{4};
% [解説 L63] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 WINDOW」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
        WINDOW=varargin{5}; NOVERLAP=varargin{6};
% [解説 L64] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 mode」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
        mode='local';   
end

% [解説 L67] 信号またはベクトルの長さを取得します。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
T=length(s_target);

% [解説 L69] 値 mode に応じて、入力形式・音源モデル・評価モードなどの処理経路を切り替えます。
switch mode        
% [解説 L70] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 'global'
% [解説 L71] 値 isempty(e_noise) に応じて、入力形式・音源モデル・評価モードなどの処理経路を切り替えます。
        switch isempty(e_noise)
% [解説 L72] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
            case 1
                % Computation of the energy ratios
% [解説 L74] この関数の複数の計算結果を受け取り、SDR、SIR、SAR にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                [SDR,SIR,SAR]=bss_energy_ratios(s_target,e_interf,e_artif);
% [解説 L75] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
                varargout{1}=10*log10(SDR); varargout{2}=10*log10(SIR); varargout{3}=10*log10(SAR);
% [解説 L76] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
            case 0
                % Computation of the energy ratios
% [解説 L78] この関数の複数の計算結果を受け取り、SDR、SIR、SNR、SAR にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                [SDR,SIR,SNR,SAR]=bss_energy_ratios(s_target,e_interf,e_noise,e_artif);
% [解説 L79] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
                varargout{1}=10*log10(SDR); varargout{2}=10*log10(SIR);
% [解説 L80] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
                varargout{3}=10*log10(SNR); varargout{4}=10*log10(SAR);                
        end
        
% [解説 L83] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 'local'
        
% [解説 L85] 信号またはベクトルの長さを取得します。結果を「各周波数で観測を音源へ分ける分離行列」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
        W=length(WINDOW); % Length of window
% [解説 L86] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 n_frames」へ代入します。
% [意図] 時間ずれまたはフィルタ効果を含む参照部分空間を作り、許容される成分を推定信号から射影するためです。
        n_frames = fix((T-NOVERLAP)/(W-NOVERLAP)); % Number of frames
        
% [解説 L71] 値 isempty(e_noise) に応じて、入力形式・音源モデル・評価モードなどの処理経路を切り替えます。
        switch isempty(e_noise)
% [解説 L89] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
            case 1
% [解説 L90] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 F_s_target」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
                F_s_target=bss_make_frames(s_target,WINDOW,NOVERLAP);
% [解説 L91] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 F_e_interf」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
                F_e_interf=bss_make_frames(e_interf,WINDOW,NOVERLAP);
% [解説 L92] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 F_e_artif」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
                F_e_artif=bss_make_frames(e_artif,WINDOW,NOVERLAP);
% [解説 L93] この関数の複数の計算結果を受け取り、SDR、SIR、SAR にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                [SDR,SIR,SAR]=bss_energy_ratios(F_s_target,F_e_interf,F_e_artif);
% [解説 L94] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
                varargout{1}=10*log10(SDR); varargout{2}=10*log10(SIR); varargout{3}=10*log10(SAR);
% [解説 L95] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
            case 0
% [解説 L96] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 F_s_target」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
                F_s_target=bss_make_frames(s_target,WINDOW,NOVERLAP);
% [解説 L97] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 F_e_interf」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
                F_e_interf=bss_make_frames(e_interf,WINDOW,NOVERLAP);
% [解説 L98] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 F_e_noise」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
                F_e_noise=bss_make_frames(e_noise,WINDOW,NOVERLAP);
% [解説 L99] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 F_e_artif」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
                F_e_artif=bss_make_frames(e_artif,WINDOW,NOVERLAP);
% [解説 L100] この関数の複数の計算結果を受け取り、SDR、SIR、SNR、SAR にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                [SDR,SIR,SNR,SAR]=bss_energy_ratios(F_s_target,F_e_interf,F_e_noise,F_e_artif);
% [解説 L101] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
                varargout{1}=10*log10(SDR); varargout{2}=10*log10(SIR);
% [解説 L102] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
                varargout{3}=10*log10(SNR); varargout{4}=10*log10(SAR);
        end        
% [解説 L104] この行で、BSS評価の成分分解または可変出力に必要な処理を行います。
% [意図] 推定信号を正解成分と誤差成分へ分け、分離性能を定量評価するためです。
end %mode
