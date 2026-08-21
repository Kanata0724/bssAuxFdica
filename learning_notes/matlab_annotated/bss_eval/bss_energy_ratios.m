% ======================================================================
% 学習用・日本語解説付きコピー
% 元ファイル: bss_eval/bss_energy_ratios.m
% 注意: 元コードは変更していません。空行と元コメントには追加解説を付けていません。
% 概要: 正解成分・干渉成分・アーティファクト成分のエネルギー比から評価指標を計算します。
% ======================================================================

% [解説 L1] この関数は入力を受け取り、varargout を計算して返します。
% [意図] ファイル全体では「正解成分・干渉成分・アーティファクト成分のエネルギー比から評価指標を計算します。」という役割を担当します。
function [varargout]=bss_energy_ratios(F_s_target,F_e_interf,varargin)

% compute energy ratios corresponding to SDR/SIR/SNR/SAR given a decomposition of an estimated source into target/interference/noise/artifacts over frames.
%
% Usage:
%
%    [SDR,SIR,SAR]    =bss_energy_ratios(F_s_target,F_e_interf,F_e_artif)
%    [SDR,SIR,SNR,SAR]=bss_energy_ratios(F_s_target,F_e_interf,F_e_noise,F_e_artif)
%
% Input:
%   - F_s_target: n_frames x T matrix containing the frames of the target source contribution,
%   - F_e_interf: n_frames x T matrix containing the frames of the interferences contribution,
%   - F_e_noise: n_frames x T matrix containing the frames of the noise contribution (if any),
%   - F_e_artif: n_frames x T matrix containing the frames of the artifacts contribution.
%
% Ouput:
%   - SDR: n_frames x 1 vector contaning the Source to Distortion Ratios per frame,
%   - SIR: n_frames x 1 vector contaning the Source to Interferences Ratios per frame,
%   - SNR: n_frames x 1 vector contaning the Signal to Noise Ratios (if noise) per frame,
%   - SAR: n_frames x 1 vector contaning the Signal to Artifacts Ratios per frame.
%
% Developers:  - Cedric Fevotte (cf269@cam.ac.uk) - Emmanuel Vincent
% (vincent@ircam.fr) - Remi Gribonval (remi.gribonval@irisa.fr)

% [解説 L25] 値 nargin に応じて、入力形式・音源モデル・評価モードなどの処理経路を切り替えます。
switch nargin
% [解説 L26] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 3
% [解説 L27] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 F_e_artif」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
        F_e_artif=varargin{1};
        % SDR
% [解説 L29] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 F_e_total」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
        F_e_total=F_e_interf+F_e_artif;
% [解説 L30] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
        varargout{1}= sum(F_s_target.^2,2)./sum(F_e_total.^2,2);
        % SIR
% [解説 L32] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
        varargout{2}=sum(F_s_target.^2,2)./sum(F_e_interf.^2,2);
        % SAR
% [解説 L34] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
        varargout{3}=sum((F_s_target+F_e_interf).^2,2)./sum(F_e_artif.^2,2);        
        
% [解説 L36] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 4        
% [解説 L37] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 F_e_noise」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
        F_e_noise=varargin{1};
% [解説 L38] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 F_e_artif」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
        F_e_artif=varargin{2};
        % SDR
% [解説 L40] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 F_e_total」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
        F_e_total=F_e_interf+F_e_noise+F_e_artif;
% [解説 L41] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
        varargout{1}=sum(F_s_target.^2,2)./sum(F_e_total.^2,2);
        % SIR
% [解説 L43] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
        varargout{2}=sum(F_s_target.^2,2)./sum(F_e_interf.^2,2);
        % SNR
% [解説 L45] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
        varargout{3}=sum((F_s_target+F_e_interf).^2,2)./sum(F_e_noise.^2,2);
        % SAR
% [解説 L47] 計算済みの成分または評価値を、要求された可変出力位置へ格納します。
% [意図] 呼び出し側が要求した評価値だけを、MATLABの可変個数出力として正しい順序で返すためです。
        varargout{4}=sum((F_s_target+F_e_interf+F_e_noise).^2,2)./sum(F_e_artif.^2,2);        
end
