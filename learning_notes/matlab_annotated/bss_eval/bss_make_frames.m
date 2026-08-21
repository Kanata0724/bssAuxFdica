% ======================================================================
% 学習用・日本語解説付きコピー
% 元ファイル: bss_eval/bss_make_frames.m
% 注意: 元コードは変更していません。空行と元コメントには追加解説を付けていません。
% 概要: BSS Evalで使うフレームまたは遅延信号の行列を組み立てます。
% ======================================================================

% [解説 L1] この関数は入力を受け取り、F_S frames_index を計算して返します。
% [意図] ファイル全体では「BSS Evalで使うフレームまたは遅延信号の行列を組み立てます。」という役割を担当します。
function [F_S frames_index] = bss_make_frames(S,WINDOW,NOVERLAP)

% decompose some signal(s) into frames
%
% Usage: [F_S frames_index] = bss_make_frames(S,WINDOW,NOVERLAP)
%
% Input:
%   - S: matrix of size n x T (with T>n),
%   - WINDOW: 1 x W window
%   - NOVERLAP: number of samples overlap
%
% Output:
%   - F_S: 
%       * if n=1, F_S is a n_frames x W matrix containing the frames (of length W) in
%       rows,
%       * if n>1, F_S is a n_frames x W x n tensor containing the frames
%       decomposition of each row of S.
%   - frames_index: index of the beginning of each frame in the rows of S
%
% Developers:  - Cedric Fevotte (cf269@cam.ac.uk) - Emmanuel Vincent
% (vincent@ircam.fr) - Remi Gribonval (remi.gribonval@irisa.fr)

% [解説 L23] この関数の複数の計算結果を受け取り、n、T にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[n,T]=size(S);

% [解説 L25] 条件「n>T」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if n>T
% [解説 L26] 実験の進行状況または計算結果を表示し、処理がどこまで進んだか確認できるようにします。
    disp('Wrong dimensions: must have T>n.')
% [解説 L27] 必要な結果が確定したため、残りの処理を行わず呼び出し元へ戻ります。
    return;
end

%%% Default values %%%
% [解説 L31] 信号またはベクトルの長さを取得します。結果を「各周波数で観測を音源へ分ける分離行列」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
W=length(WINDOW); % Length of window

% [解説 L33] 条件「T < W」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if T < W
% [解説 L34] 実験の進行状況または計算結果を表示し、処理がどこまで進んだか確認できるようにします。
    disp('Please choose a window smaller than the signals.')
% [解説 L35] 必要な結果が確定したため、残りの処理を行わず呼び出し元へ戻ります。
    return;
end

% [解説 L38] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 n_frames」へ代入します。
% [意図] 時間ずれまたはフィルタ効果を含む参照部分空間を作り、許容される成分を推定信号から射影するためです。
n_frames = fix((T-NOVERLAP)/(W-NOVERLAP)); % Number of frames
% If needed the very end of the signal is removed.

% [解説 L41] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 frames_index」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
frames_index = 1 + (0:(n_frames-1))*(W-NOVERLAP); % Index of beginnings of frames

% [解説 L43] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 F_S」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
F_S=zeros(n_frames,W,n); % If n=1, F_S is a 2-D tensor (matrix)

% [解説 L45] 処理途中で使う変数 i を 1:n の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 全周波数ビンへ同じ処理を適用し、帯域の一部だけが未処理になるのを防ぐためです。
for i=1:n
% [解説 L46] 処理途中で使う変数 k を 1:n_frames の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
    for k=1:n_frames
% [解説 L47] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 F_S(k,:,i)」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
        F_S(k,:,i)=(S(i,frames_index(k)+(0:(W-1))).*WINDOW);
    end
end
