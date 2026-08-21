% ======================================================================
% 学習用・日本語解説付きコピー
% 元ファイル: bss_eval/bss_eval_images.m
% 注意: 元コードは変更していません。空行と元コメントには追加解説を付けていません。
% 概要: マルチチャネルの参照音像と推定音像を対応付け、分離性能を評価します。
% ======================================================================

% [解説 L1] この関数は入力を受け取り、SDR,ISR,SIR,SAR,perm を計算して返します。
% [意図] ファイル全体では「マルチチャネルの参照音像と推定音像を対応付け、分離性能を評価します。」という役割を担当します。
function [SDR,ISR,SIR,SAR,perm]=bss_eval_images(ie,i)

% BSS_EVAL_IMAGES Ordering and measurement of the separation quality for
% estimated source spatial image signals in terms of true source, spatial
% (or filtering) distortion, interference and artifacts.
%
% [SDR,ISR,SIR,SAR,perm]=bss_eval_images(ie,i)
%
% Inputs:
% ie: nsrc x nsampl x nchan matrix containing estimated source images
% i: nsrc x nsampl x nchan matrix containing true source images
%
% Outputs:
% SDR: nsrc x 1 vector of Signal to Distortion Ratios
% ISR: nsrc x 1 vector of source Image to Spatial distortion Ratios
% SIR: nsrc x 1 vector of Source to Interference Ratios
% SAR: nsrc x 1 vector of Sources to Artifacts Ratios
% perm: nsrc x 1 vector containing the best ordering of estimated source
% images in the mean SIR sense (estimated source image number perm(j)
% corresponds to true source image number j)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2007-2008 Emmanuel Vincent
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
% If you find it useful, please cite the following reference:
% Emmanuel Vincent, Hiroshi Sawada, Pau Bofill, Shoji Makino and Justinian
% P. Rosca, "First stereo audio source separation evaluation campaign:
% data, algorithms and results," In Proc. Int. Conf. on Independent
% Component Analysis and Blind Source Separation (ICA), pp. 552-559, 2007.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Errors %%%
% [解説 L35] 条件「nargin<2, error('Not enough input arguments.'); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if nargin<2, error('Not enough input arguments.'); end
% [解説 L36] この関数の複数の計算結果を受け取り、nsrc、nsampl、nchan にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[nsrc,nsampl,nchan]=size(ie);
% [解説 L37] この関数の複数の計算結果を受け取り、nsrc2、nsampl2、nchan2 にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[nsrc2,nsampl2,nchan2]=size(i);
% [解説 L38] 条件「nsrc2~=nsrc, error('The number of estimated source images and reference source images must be equal.'); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if nsrc2~=nsrc, error('The number of estimated source images and reference source images must be equal.'); end
% [解説 L39] 条件「nsampl2~=nsampl, error('The estimated source images and reference source images must have the same duration.'); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if nsampl2~=nsampl, error('The estimated source images and reference source images must have the same duration.'); end
% [解説 L40] 条件「nchan2~=nchan, error('The estimated source images and reference source images must have the same number of channels.'); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if nchan2~=nchan, error('The estimated source images and reference source images must have the same number of channels.'); end

%%% Performance criteria %%%
% Computation of the criteria for all possible pair matches
% [解説 L44] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「最適化または分離性能を表す評価値」へ代入します。
% [意図] 分離品質を目的別の数値で比較できるようにするためです。
SDR=zeros(nsrc,nsrc);
% [解説 L45] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「最適化または分離性能を表す評価値」へ代入します。
% [意図] 分離品質を目的別の数値で比較できるようにするためです。
ISR=zeros(nsrc,nsrc);
% [解説 L46] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「最適化または分離性能を表す評価値」へ代入します。
% [意図] 分離品質を目的別の数値で比較できるようにするためです。
SIR=zeros(nsrc,nsrc);
% [解説 L47] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「最適化または分離性能を表す評価値」へ代入します。
% [意図] 分離品質を目的別の数値で比較できるようにするためです。
SAR=zeros(nsrc,nsrc);
% [解説 L48] 処理途中で使う変数 jest を 1:nsrc, の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
for jest=1:nsrc,
% [解説 L49] 処理途中で使う変数 jtrue を 1:nsrc, の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
    for jtrue=1:nsrc,
% [解説 L50] この関数の複数の計算結果を受け取り、s_true、e_spat、e_interf、e_artif にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
        [s_true,e_spat,e_interf,e_artif]=bss_decomp_mtifilt(reshape(ie(jest,:,:),nsampl,nchan).',i,jtrue,512);
% [解説 L51] この関数の複数の計算結果を受け取り、SDR(jest,jtrue)、ISR(jest,jtrue)、SIR(jest,jtrue)、SAR(jest,jtrue) にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
        [SDR(jest,jtrue),ISR(jest,jtrue),SIR(jest,jtrue),SAR(jest,jtrue)]=bss_image_crit(s_true,e_spat,e_interf,e_artif);
    end
end
% Selection of the best ordering
% [解説 L55] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「音源番号の対応を揃える置換情報」へ代入します。
% [意図] 周波数ごとに入れ替わり得る音源番号を統一し、時間波形として復元できるようにするためです。
perm=perms(1:nsrc);
% [解説 L56] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 nperm」へ代入します。
% [意図] 音源特徴の類似度と候補順を保持し、周波数間の音源番号を一貫した順序へ揃えるためです。
nperm=size(perm,1);
% [解説 L57] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 meanSIR」へ代入します。
% [意図] 分離品質をエネルギー比として数値化し、入力条件や手法間で公平に比較するためです。
meanSIR=zeros(nperm,1);
% [解説 L58] 処理途中で使う変数 p を 1:nperm, の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
for p=1:nperm,
% [解説 L59] 指定した軸で平均し、ばらつきの少ない代表量を作ります。結果を「処理途中で使う変数 meanSIR(p)」へ代入します。
% [意図] 分離品質をエネルギー比として数値化し、入力条件や手法間で公平に比較するためです。
    meanSIR(p)=mean(SIR((0:nsrc-1)*nsrc+perm(p,:)));
end
% [解説 L61] この関数の複数の計算結果を受け取り、meanSIR、popt にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[meanSIR,popt]=max(meanSIR);
% [解説 L62] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「音源番号の対応を揃える置換情報」へ代入します。
% [意図] 周波数ごとに入れ替わり得る音源番号を統一し、時間波形として復元できるようにするためです。
perm=perm(popt,:).';
% [解説 L63] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「最適化または分離性能を表す評価値」へ代入します。
% [意図] 分離品質を目的別の数値で比較できるようにするためです。
SDR=SDR((0:nsrc-1).'*nsrc+perm);
% [解説 L64] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「最適化または分離性能を表す評価値」へ代入します。
% [意図] 分離品質を目的別の数値で比較できるようにするためです。
ISR=ISR((0:nsrc-1).'*nsrc+perm);
% [解説 L65] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「最適化または分離性能を表す評価値」へ代入します。
% [意図] 分離品質を目的別の数値で比較できるようにするためです。
SIR=SIR((0:nsrc-1).'*nsrc+perm);
% [解説 L66] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「最適化または分離性能を表す評価値」へ代入します。
% [意図] 分離品質を目的別の数値で比較できるようにするためです。
SAR=SAR((0:nsrc-1).'*nsrc+perm);

% [解説 L68] 必要な結果が確定したため、残りの処理を行わず呼び出し元へ戻ります。
return;



% [解説 L72] この関数は入力を受け取り、s_true,e_spat,e_interf,e_artif を計算して返します。
% [意図] ファイル全体では「マルチチャネルの参照音像と推定音像を対応付け、分離性能を評価します。」という役割を担当します。
function [s_true,e_spat,e_interf,e_artif]=bss_decomp_mtifilt(se,s,j,flen)

% BSS_DECOMP_MTIFILT Decomposition of an estimated source image into four
% components representing respectively the true source image, spatial (or
% filtering) distortion, interference and artifacts, derived from the true
% source images using multichannel time-invariant filters.
%
% [s_true,e_spat,e_interf,e_artif]=bss_decomp_mtifilt(se,s,j,flen)
%
% Inputs:
% se: nchan x nsampl matrix containing the estimated source image (one row per channel)
% s: nsrc x nsampl x nchan matrix containing the true source images
% j: source index corresponding to the estimated source image in s
% flen: length of the multichannel time-invariant filters in samples
%
% Outputs:
% s_true: nchan x nsampl matrix containing the true source image (one row per channel)
% e_spat: nchan x nsampl matrix containing the spatial (or filtering) distortion component
% e_interf: nchan x nsampl matrix containing the interference component
% e_artif: nchan x nsampl matrix containing the artifacts component

%%% Errors %%%
% [解説 L94] 条件「nargin<4, error('Not enough input arguments.'); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if nargin<4, error('Not enough input arguments.'); end
% [解説 L95] この関数の複数の計算結果を受け取り、nchan2、nsampl2 にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[nchan2,nsampl2]=size(se);
% [解説 L96] この関数の複数の計算結果を受け取り、nsrc、nsampl、nchan にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[nsrc,nsampl,nchan]=size(s);
% [解説 L97] 条件「nchan2~=nchan, error('The number of channels of the true source images and the estimated source image must be equal.'); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if nchan2~=nchan, error('The number of channels of the true source images and the estimated source image must be equal.'); end
% [解説 L98] 条件「nsampl2~=nsampl, error('The duration of the true source images and the estimated source image must be equal.'); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if nsampl2~=nsampl, error('The duration of the true source images and the estimated source image must be equal.'); end

%%% Decomposition %%%
% True source image
% [解説 L102] 要素数を保ったままshapeを組み替えます。結果を「処理途中で使う変数 s_true」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
s_true=[reshape(s(j,:,:),nsampl,nchan).',zeros(nchan,flen-1)];
% Spatial (or filtering) distortion
% [解説 L104] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 e_spat」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
e_spat=project(se,s(j,:,:),flen)-s_true;
% Interference
% [解説 L106] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 e_interf」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
e_interf=project(se,s,flen)-s_true-e_spat;
% Artifacts
% [解説 L108] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 e_artif」へ代入します。
% [意図] 推定誤差を正解・干渉・雑音・アーティファクトへ分け、性能低下の原因を区別するためです。
e_artif=[se,zeros(nchan,flen-1)]-s_true-e_spat-e_interf;

% [解説 L110] 必要な結果が確定したため、残りの処理を行わず呼び出し元へ戻ります。
return;



% [解説 L114] この関数は入力を受け取り、sproj を計算して返します。
% [意図] ファイル全体では「マルチチャネルの参照音像と推定音像を対応付け、分離性能を評価します。」という役割を担当します。
function sproj=project(se,s,flen)

% SPROJ Least-squares projection of each channel of se on the subspace
% spanned by delayed versions of the channels of s, with delays between 0
% and flen-1

% [解説 L120] この関数の複数の計算結果を受け取り、nsrc、nsampl、nchan にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[nsrc,nsampl,nchan]=size(s);
% [解説 L121] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「正解音源または正解音像の時間信号」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
s=reshape(permute(s,[3 1 2]),nchan*nsrc,nsampl);

%%% Computing coefficients of least squares problem via FFT %%%
% Zero padding and FFT of input data
% [解説 L125] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「正解音源または正解音像の時間信号」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
s=[s,zeros(nchan*nsrc,flen-1)];
% [解説 L126] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「FDICAが推定した分離時間信号」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
se=[se,zeros(nchan,flen-1)];
% [解説 L127] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 fftlen」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
fftlen=2^nextpow2(nsampl+flen-1);
% [解説 L128] 時間方向の信号を周波数成分へ変換します。結果を「処理途中で使う変数 sf」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
sf=fft(s,fftlen,2);
% [解説 L129] 時間方向の信号を周波数成分へ変換します。結果を「処理途中で使う変数 sef」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
sef=fft(se,fftlen,2);
% Inner products between delayed versions of s
% [解説 L131] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 G」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
G=zeros(nchan*nsrc*flen);
% [解説 L132] 処理途中で使う変数 k1 を 0:nchan*nsrc-1, の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
for k1=0:nchan*nsrc-1,
% [解説 L133] 処理途中で使う変数 k2 を 0:k1, の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
    for k2=0:k1,
% [解説 L134] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 ssf」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
        ssf=sf(k1+1,:).*conj(sf(k2+1,:));
% [解説 L135] 周波数成分を時間方向の信号へ戻します。結果を「処理途中で使う変数 ssf」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
        ssf=real(ifft(ssf));
% [解説 L136] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 ss」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
        ss=toeplitz(ssf([1 fftlen:-1:fftlen-flen+2]),ssf(1:flen));
% [解説 L137] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 G(k1*flen+1:k1*flen+flen,k2*flen+1:k2*flen+flen)」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
        G(k1*flen+1:k1*flen+flen,k2*flen+1:k2*flen+flen)=ss;
% [解説 L138] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 G(k2*flen+1:k2*flen+flen,k1*flen+1:k1*flen+flen)」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
        G(k2*flen+1:k2*flen+flen,k1*flen+1:k1*flen+flen)=ss.';
    end
end
% Inner products between se and delayed versions of s
% [解説 L142] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 D」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
D=zeros(nchan*nsrc*flen,nchan);
% [解説 L143] 処理途中で使う変数 k を 0:nchan*nsrc-1, の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
for k=0:nchan*nsrc-1,
% [解説 L144] 処理途中で使う変数 i を 1:nchan, の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 全周波数ビンへ同じ処理を適用し、帯域の一部だけが未処理になるのを防ぐためです。
    for i=1:nchan,
% [解説 L145] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 ssef」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
        ssef=sf(k+1,:).*conj(sef(i,:));
% [解説 L146] 周波数成分を時間方向の信号へ戻します。結果を「処理途中で使う変数 ssef」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
        ssef=real(ifft(ssef,[],2));
% [解説 L147] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 D(k*flen+1:k*flen+flen,i)」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
        D(k*flen+1:k*flen+flen,i)=ssef(:,[1 fftlen:-1:fftlen-flen+2]).';
    end
end

%%% Computing projection %%%
% Distortion filters
% [解説 L153] 逆行列を明示せず線形方程式を解き、数値的に安定な係数を求めます。結果を「処理途中で使う変数 C」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
C=G\D;
% [解説 L154] 要素数を保ったままshapeを組み替えます。結果を「処理途中で使う変数 C」へ代入します。
% [意図] 正解との対応付けと誤差成分分解を進め、最終的なSDR・SIR・SARへ結び付けるためです。
C=reshape(C,flen,nchan*nsrc,nchan);
% Filtering
% [解説 L156] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 sproj」へ代入します。
% [意図] 時間ずれまたはフィルタ効果を含む参照部分空間を作り、許容される成分を推定信号から射影するためです。
sproj=zeros(nchan,nsampl+flen-1);
% [解説 L157] 処理途中で使う変数 k を 1:nchan*nsrc, の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 対象配列の全要素へ同じ規則を適用し、完全な結果を組み立てるためです。
for k=1:nchan*nsrc,
% [解説 L158] 処理途中で使う変数 i を 1:nchan, の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 全周波数ビンへ同じ処理を適用し、帯域の一部だけが未処理になるのを防ぐためです。
    for i=1:nchan,
% [解説 L159] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 sproj(i,:)」へ代入します。
% [意図] 時間ずれまたはフィルタ効果を含む参照部分空間を作り、許容される成分を推定信号から射影するためです。
        sproj(i,:)=sproj(i,:)+fftfilt(C(:,k,i).',s(k,:));
    end
end

% [解説 L163] 必要な結果が確定したため、残りの処理を行わず呼び出し元へ戻ります。
return;



% [解説 L167] この関数は入力を受け取り、SDR,ISR,SIR,SAR を計算して返します。
% [意図] ファイル全体では「マルチチャネルの参照音像と推定音像を対応付け、分離性能を評価します。」という役割を担当します。
function [SDR,ISR,SIR,SAR]=bss_image_crit(s_true,e_spat,e_interf,e_artif)

% BSS_IMAGE_CRIT Measurement of the separation quality for a given source
% image in terms of true source, spatial (or filtering) distortion,
% interference and artifacts.
%
% [SDR,ISR,SIR,SAR]=bss_image_crit(s_true,e_spat,e_interf,e_artif)
%
% Inputs:
% s_true: nchan x nsampl matrix containing the true source image (one row per channel)
% e_spat: nchan x nsampl matrix containing the spatial (or filtering) distortion component
% e_interf: nchan x nsampl matrix containing the interference component
% e_artif: nchan x nsampl matrix containing the artifacts component
%
% Outputs:
% SDR: Signal to Distortion Ratio
% ISR: source Image to Spatial distortion Ratio
% SIR: Source to Interference Ratio
% SAR: Sources to Artifacts Ratio

%%% Errors %%%
% [解説 L188] 条件「nargin<4, error('Not enough input arguments.'); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if nargin<4, error('Not enough input arguments.'); end
% [解説 L189] この関数の複数の計算結果を受け取り、nchant、nsamplt にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[nchant,nsamplt]=size(s_true);
% [解説 L190] この関数の複数の計算結果を受け取り、nchans、nsampls にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[nchans,nsampls]=size(e_spat);
% [解説 L191] この関数の複数の計算結果を受け取り、nchani、nsampli にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[nchani,nsampli]=size(e_interf);
% [解説 L192] この関数の複数の計算結果を受け取り、nchana、nsampla にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[nchana,nsampla]=size(e_artif);
% [解説 L193] 条件「~((nchant==nchans)&&(nchant==nchani)&&(nchant==nchana)), error('All the components must have the same number of channels.'); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if ~((nchant==nchans)&&(nchant==nchani)&&(nchant==nchana)), error('All the components must have the same number of channels.'); end
% [解説 L194] 条件「~((nsamplt==nsampls)&&(nsamplt==nsampli)&&(nsamplt==nsampla)), error('All the components must have the same duration.'); end」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if ~((nsamplt==nsampls)&&(nsamplt==nsampli)&&(nsamplt==nsampla)), error('All the components must have the same duration.'); end

%%% Energy ratios %%%
% SDR
% [解説 L198] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「最適化または分離性能を表す評価値」へ代入します。
% [意図] 分離品質を目的別の数値で比較できるようにするためです。
SDR=10*log10(sum(sum(s_true.^2))/sum(sum((e_spat+e_interf+e_artif).^2)));
% ISR
% [解説 L200] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「最適化または分離性能を表す評価値」へ代入します。
% [意図] 分離品質を目的別の数値で比較できるようにするためです。
ISR=10*log10(sum(sum(s_true.^2))/sum(sum(e_spat.^2)));
% SIR
% [解説 L202] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「最適化または分離性能を表す評価値」へ代入します。
% [意図] 分離品質を目的別の数値で比較できるようにするためです。
SIR=10*log10(sum(sum((s_true+e_spat).^2))/sum(sum(e_interf.^2)));
% SAR
% [解説 L204] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「最適化または分離性能を表す評価値」へ代入します。
% [意図] 分離品質を目的別の数値で比較できるようにするためです。
SAR=10*log10(sum(sum((s_true+e_spat+e_interf).^2))/sum(sum(e_artif.^2)));

% [解説 L206] 必要な結果が確定したため、残りの処理を行わず呼び出し元へ戻ります。
return;
