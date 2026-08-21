% ======================================================================
% 学習用・日本語解説付きコピー
% 元ファイル: bss_eval/bss_decomp_tvgain.m
% 注意: 元コードは変更していません。空行と元コメントには追加解説を付けていません。
% 概要: 推定信号を正解・干渉・アーティファクト成分へ分解します。
% ======================================================================

% [解説 L1] この関数は入力を受け取り、varargout を計算して返します。
% [意図] ファイル全体では「推定信号を正解・干渉・アーティファクト成分へ分解します。」という役割を担当します。
function varargout=bss_decomp_tvgain(varargin)

% decompose an estimated source into target/interference/noise/artefacts components, assuming the admissible distortion is a time-varying gain.
%
% Usage:
%
% [s_target,e_interf[,e_noise],e_artif]=bss_decomp_tvgain(se,index,S[,N],tvshape,tvstep)
%
% Input:
%   - se: row vector of length T containing the estimated source,
%   - index: points which component of S se has to be compared to,
%   - S: n x T matrix containing the original sources,
%   - N: m x T matrix containing the noise on the obseravtions (if any).
%   - tvshape : row vector of length V at most T containing the shape of the elementary 
%     allowed time variations of the gain
%   - tvstep  : hop size (in number of samples) between two consecutive
%     variations of the gain
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
       
% [解説 L32] 値 nargin に応じて、入力形式・音源モデル・評価モードなどの処理経路を切り替えます。
switch nargin
% [解説 L33] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 5
% [解説 L34] この関数の複数の計算結果を受け取り、varargout、varargout、varargout にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
        [varargout{1},varargout{2},varargout{3}]=bss_decomp_tvfilt(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},0);
% [解説 L35] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
    case 6
% [解説 L36] この関数の複数の計算結果を受け取り、varargout、varargout、varargout、varargout にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
        [varargout{1},varargout{2},varargout{3},varargout{4}]=bss_decomp_filt(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6},0);
% [解説 L37] 既知の選択肢に当てはまらない入力を処理し、想定外の設定を見逃さないようにします。
    otherwise
% [解説 L38] 実験の進行状況または計算結果を表示し、処理がどこまで進んだか確認できるようにします。
        disp('Wrong number of arguments.')
end

