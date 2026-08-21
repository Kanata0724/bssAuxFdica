% ======================================================================
% 学習用・日本語解説付きコピー
% 元ファイル: bss_eval/bss_proj.m
% 注意: 元コードは変更していません。空行と元コメントには追加解説を付けていません。
% 概要: 推定信号を参照信号が張る部分空間へ射影します。
% ======================================================================

% [解説 L1] この関数は入力を受け取り、PY_x coeff を計算して返します。
% [意図] ファイル全体では「推定信号を参照信号が張る部分空間へ射影します。」という役割を担当します。
function [PY_x coeff]=bss_proj(x,Y)

% compute the orthogonal projection of x on the subspace spanned by the row(s) of Y.
%
% Usage: PY_x         = proj(x,Y)
%        [PY_x coeff] = proj(x,Y)
%
% Input:
%   - x: row vector of length T,
%   - Y: vector or matrix of length T.
%
% Ouput:
%   - PY_x: row vector of length T containing the orthogonal projection of
%   x onto the range of the rows of Y.
%   - coeff : column vector with as many rows as Y containing the
%   coefficients such that PY_x = coeff.'*Y
%
% Developers:  - Cedric Fevotte (cf269@cam.ac.uk) - Emmanuel Vincent
% (vincent@ircam.fr) - Remi Gribonval (remi.gribonval@irisa.fr)

% Gram matrix of Y
% [解説 L22] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 G」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
G=Y*Y';

%same as coeff=inv(conj(G))*conj(Y*x');
% [解説 L25] 逆行列を明示せず線形方程式を解き、数値的に安定な係数を求めます。結果を「処理途中で使う変数 coeff」へ代入します。
% [意図] 時間ずれまたはフィルタ効果を含む参照部分空間を作り、許容される成分を推定信号から射影するためです。
coeff=conj(G)\conj(Y*x');
%if the Gram matrix G is not invertible then coeff=pinv(conj(G))*conj(Y*x')
%should work, but in general it is much slower than the default code

% [解説 L29] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 PY_x」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
PY_x=  coeff.'*Y;

% Same as PY_x= x*pinv(Y)*Y;
