% ======================================================================
% 学習用・日本語解説付きコピー
% 元ファイル: bss_eval/bss_tvproj.m
% 注意: 元コードは変更していません。空行と元コメントには追加解説を付けていません。
% 概要: 時間変化するフィルタを許した参照部分空間へ推定信号を射影します。
% ======================================================================

% [解説 L1] この関数は入力を受け取り、PY_x coeff を計算して返します。
% [意図] ファイル全体では「時間変化するフィルタを許した参照部分空間へ推定信号を射影します。」という役割を担当します。
function [PY_x coeff] = bss_tvproj(x,Y,tvshape,tvstep)

% compute the orthogonal projection of x onto the space of shifted windowed versions of the row(s) of Y
%
% Usage: [PY_x coeff] = tvproj(x,Y,tvshape,tvstep)
%
% Input:
%   - x: row vector of length T corresponding to the signal to be projected,
%   - Y: vector or matrix of length T with n_rows rows which windowed rows span the projection space.
%   - tvshape : row vector of length V corresponding to the window applied
%   to Y to define the projection space
%   - tvstep  : number of samples between two adjacent windows
%
% Ouput:
%   - PY_x: row vector of length T containing the orthogonal projection of
%   x onto the range of the shifted windowed versions of the row(s) of Y.
%   - coeff : matrix with n_rows rows and n_frames columns containing the coefficients 
%   of the projection 
%
% Developers:  - Cedric Fevotte (cf269@cam.ac.uk) - Emmanuel Vincent
% (vincent@ircam.fr) - Remi Gribonval (remi.gribonval@irisa.fr)

% 1. Decompose Y into frames using tvshape as a window
% [解説 L24] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「AuxFDICA更新や射影に使う共分散・相関行列」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
V        = size(tvshape,2); % the size of the window
% [解説 L25] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 NOVERLAP」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
NOVERLAP = V-tvstep;        % convert the hop size into a number of overlapping samples

% [解説 L27] この関数の複数の計算結果を受け取り、Y_frames、frames_index にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[Y_frames frames_index] = bss_make_frames(Y,tvshape,NOVERLAP); % Y_frames is a 3-D array
% [解説 L28] この関数の複数の計算結果を受け取り、n_frames、V1、n_rows にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[n_frames V1 n_rows]= size(Y_frames);

% 2. Compute the inner products between x and the frames of the row(s) of Y
% [解説 L31] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 ip」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
ip = zeros(n_frames,n_rows);
% [解説 L32] 現在処理している周波数ビン番号 を 1:n_frames % loop on frames の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
for f=1:n_frames % loop on frames
% [解説 L33] 要素数を保ったままshapeを組み替えます。結果を「処理途中で使う変数 ip(f,:)」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
    ip(f,:) = x(frames_index(f)+(0:(V-1)) ) * reshape(Y_frames(f,:,:),V,n_rows); % columns of ip correspond to the same frame number
end
% [解説 L35] 要素数を保ったままshapeを組み替えます。結果を「処理途中で使う変数 ip」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
ip = reshape(ip,n_frames*n_rows,1); % now we want a single column vector to apply an inverse matrix to it

% 3. Compute the Gram matrix, which is square of size (n_frames x n_rows) x
% (n_frames x n_rows) 
% [解説 L39] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 Gram」へ代入します。
% [意図] 時間ずれまたはフィルタ効果を含む参照部分空間を作り、許容される成分を推定信号から射影するためです。
Gram = zeros(n_frames*n_rows,n_frames*n_rows);
% [解説 L40] 現在処理している周波数ビン番号 を 1:n_frames の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
for f=1:n_frames 
% [解説 L41] 処理途中で使う変数 f1 を 1:n_frames の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
    for f1 = 1:n_frames
        % locate the range of the intersection between frames
% [解説 L43] 候補の最大値を取り、ピークや最良候補を求めます。結果を「処理途中で使う変数 first」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
        first = max(frames_index(f),frames_index(f1));
% [解説 L44] 候補の最小値を取り、最小コストの候補を求めます。結果を「処理途中で使う変数 last」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
        last  = min(frames_index(f),frames_index(f1))+V-1;
        % if the intersection is non empty, fill it
% [解説 L46] 条件「(first <= last)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
        if (first <= last)
% [解説 L47] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 trange」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
            trange = (first:last)-frames_index(f)+1;
% [解説 L48] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 trange1」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
            trange1= (first:last)-frames_index(f1)+1;
% [解説 L49] 要素数を保ったままshapeを組み替えます。結果を「処理途中で使う変数 Gram((f-1)*n_rows+(1:n_rows),(f1-1)*n_rows+(1:n_rows))」へ代入します。
% [意図] 時間ずれまたはフィルタ効果を含む参照部分空間を作り、許容される成分を推定信号から射影するためです。
            Gram((f-1)*n_rows+(1:n_rows),(f1-1)*n_rows+(1:n_rows)) = reshape(Y_frames(f,trange,:),length(trange),n_rows)'*reshape(Y_frames(f1,trange1,:),length(trange1),n_rows); % shall we reshape Y_frames ????
        end
    end
end

% 4. Apply the inverse of the Gram matrix to ip to get coeff
% [解説 L55] 逆行列を明示せず線形方程式を解き、数値的に安定な係数を求めます。結果を「処理途中で使う変数 coeff」へ代入します。
% [意図] 時間ずれまたはフィルタ効果を含む参照部分空間を作り、許容される成分を推定信号から射影するためです。
coeff = Gram\ip;
% [解説 L56] 要素数を保ったままshapeを組み替えます。結果を「処理途中で使う変数 coeff」へ代入します。
% [意図] 時間ずれまたはフィルタ効果を含む参照部分空間を作り、許容される成分を推定信号から射影するためです。
coeff = reshape(coeff,n_frames,n_rows);
% 4.  Reconstruct using coefficients alpha
% [解説 L58] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 PY_x」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
PY_x = zeros(size(x));

% [解説 L60] 現在処理している周波数ビン番号 を 1:n_frames の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
for f=1:n_frames
% [解説 L61] 要素数を保ったままshapeを組み替えます。結果を「処理途中で使う変数 PY_x(frames_index(f)+(0:(V-1)))」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
    PY_x(frames_index(f)+(0:(V-1))) = PY_x(frames_index(f)+(0:(V-1))) + coeff(f,:)*reshape(Y_frames(f,:,:),V,n_rows)';
end

