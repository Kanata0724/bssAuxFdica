% ======================================================================
% 学習用・日本語解説付きコピー
% 元ファイル: DGTtool.m
% 注意: 元コードは変更していません。空行と元コメントには追加解説を付けていません。
% 概要: 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
% ======================================================================

% [解説 L1] DGT/STFTの設定と内部状態を保持するDGTtoolクラスを定義します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
classdef DGTtool < handle
    %DGTtool: A simple and user-friendly tool for computing STFT/DGT.
    %   (MATLAB 2020b or later and Signal Processing Toolbox are required.)
    %
    %
    %   --- Quick start ---
    %
    %   DGTtool provides an easy way to compute STFT/DGT. 
    %    - 1st step: Construct a DGTtool object.
    %         F = DGTtool
    %    - 2nd step: Compute spectrogram X from a time-domain signal x.
    %         X = F(x);
    %   That's all for computing STFT/DGT.
    %
    %   The following command reconstructs the time-domain signal.
    %         x = F.pinv(X);
    %   Note: The reconstructed signal may be longer than the original one.
    %
    %
    %   --- Parameter settings ---
    %
    %   Parameters of STFT/DGT can be set as follows:
    %      F = DGTtool('windowShift',100,'windowLength',1000,'FFTnum',500)
    %
    %   The acceptable parameters are
    %      'windowShift'  (positive integer)
    %      'windowLength' (positive integer)
    %      'FFTnum'       (positive integer)
    %      'windowName'   (char or string)
    %      'windowVector' (column vector)
    %
    %   List of window names is given by typing DGTtool.windowList at the command line.
    %   Partial matching of leading characters is supported (case-insensitive).
    %
    %   The order of parameters can be altered as follows:
    %      F = DGTtool('windowName','b','windowLength',1000,'windowShift',100)
    %
    %   Note: MATLAB 2021a or later allows the Name=Value syntax:
    %      F = DGTtool(windowName='b', windowLength=1000, windowShift=100)
    %
    %
    %   --- Methods ---
    %
    %   The following methods can be called as F.methodName(inputVars).
    %
    %   DGTtool Methods:
    %   Forward Transforms
    %      subsref  - Implemented for shortcut notation of DGT (F(x))
    %      DGT      - Compute STFT/DGT (F(x) is shortcut of F.DGT(x))
    %      reassign - Compute sparse (reassigned) spectrogram
    %
    %   Inverse transforms
    %      H    - Inverse STFT/DGT (complex conjugate transpose of F)
    %      pinv - Inverse STFT/DGT with perfect reconstruction (pseudo-inverse of F)
    %
    %   Plot functions
    %      plot         - Draw spectrogram
    %      plotPhase    - Visualize phase spectrogram
    %      plotReassign - Draw sparse (reassigned) spectrogram
    %
    %   Window utilities
    %      setWindow       - Change window
    %      makeWindowTight - Compute canonical tight window
    %      plotWin         - Draw windows
    %
    %   Phase manipulation
    %      makeZeroPhase     - Remove linear phase component of window
    %      undoMakeZeroPhase - Cancel the effect of makeZeroPhase
    %      changeDGTdef      - Convert DGT definition
    %      undoChangeDGTdef  - Cancel the effect of changeDGTdef
    %
    %
    %   --- Static methods ---
    %
    %   Static methods can be used without constructing a DGTtool object.
    %   To use them, type DGTtool.methodName with input/output arguments.
    %
    %   DGTtool Methods:
    %   Window utilities
    %      windowList              - Returns acceptable window names
    %      getWindow               - Compute window
    %      computeCanonicalDual    - Compute canonical dual window
    %      computeCanonicalTight   - Compute canonical tight window
    %      computeNumericalDiffWin - Compute numerical differential
    %      isdual                  - Check whether a pair of windows is dual
    %
    %   Zero-padding
    %      zeroPad               - Add zero at the end of signals
    %      extendSignalByZeroPad - Add zero at the end of signals
    %      zeroPadForFactorDGT   - Add zero at the end of signals
    %      sigLenForFactorDGT    - Compute required signal length
    %
    %
    %   See the associated demo file (demo.m) for explanation and usage.
    %   See also DGTtool, DGT, pinv, plot
    
    %   Author: Kohei Yatabe (2021)
    
    % [Memo]
    % 2020b is required because of pagemtimes.
    % 2019b might be enough for performing STFT/DGT (depending on settings).
    % 2019b is required because of arguments.
    % Signal Processing Toolbox is required because of buffer.
    % Signal Processing Toolbox is required for Slepian and Chebyshev windows.
    
% [解説 L106] 窓、シフト量、FFT数、変換用添字など、このクラスが共有する状態を宣言します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
    properties (Dependent)
% [解説 L107] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        redundancy % Ratio of number of time-frequency bins to signal length.
    end
    
% [解説 L110] 窓、シフト量、FFT数、変換用添字など、このクラスが共有する状態を宣言します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
    properties
% [解説 L111] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        shift   % Amount of shift of window
% [解説 L112] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        FFTnum  % Number of frequency bins in time-frequency domain
% [解説 L113] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        win     % Column vector of window
% [解説 L114] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        dualWin % Column vector of dual window of win
    end
    
% [解説 L117] 窓、シフト量、FFT数、変換用添字など、このクラスが共有する状態を宣言します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
    properties (SetAccess = private)
% [解説 L118] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        diffWin             % Differential of win
% [解説 L119] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        sigLen              % Length of signal that is inputted last time
% [解説 L120] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 isDual」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
        isDual = false      % True when dualWin is dual of win
% [解説 L121] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 isCanonical」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
        isCanonical = false % True when dualWin is canonical dual of win
    end
    
% [解説 L124] 窓、シフト量、FFT数、変換用添字など、このクラスが共有する状態を宣言します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
    properties (Hidden = true)
% [解説 L125] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        OLAindex
% [解説 L126] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        factorIdx
% [解説 L127] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        defConverter
% [解説 L128] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        zeroPhaseConverter
% [解説 L129] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 isNotCompDual」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
        isNotCompDual = true
% [解説 L130] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 isNotWinCalcInDGT」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
        isNotWinCalcInDGT = true
    end
    
% [解説 L133] DGT、逆変換、描画、窓生成、設定更新を行うメソッド群を定義します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
    methods % main
        
% [解説 L135] この関数は入力を受け取り、obj を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function obj = DGTtool(options)
            %Constructor: Create a DGTtool object with given parameters.
            %   Parameters of STFT/DGT can be set through Name-Value pairs.
            %
            %   Usage:
            %      F = DGTtool
            %      F = DGTtool(Name,Value)
            %      F = DGTtool(Name1,Value1,Name2,Value2,...)
            %
            %   Name of options:
            %      'windowShift'  (positive integer, default = 256)
            %      'windowLength' (positive integer, default = 2048)
            %      'FFTnum'       (positive integer, default = length(win))
            %      'windowName'   (char or string,   default = '4termC5Nuttall')
            %      'windowVector' (column vector)
            %
            %   Acceptable window name:
            %      'Hann'
            %      'Blackman'
            %      '3termC1Nuttall'
            %      '3termC3Nuttall'
            %      '4termC1Nuttall'
            %      '4termC3Nuttall'
            %      '4termC5Nuttall'
            %      'Gauss'
            %      'Slepian'
            %      'Chebyshev'
            %
            %   Note: If windowVector is given, windowLength and windowName are ignored.
            %
            %   See also windowList, setWindow, plotWin
% [解説 L166] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L167] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                options.windowShift  (1,1) {mustBePositive,mustBeInteger} = 256
% [解説 L168] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                options.windowLength (1,1) {mustBePositive,mustBeInteger} = 2048
% [解説 L169] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                options.FFTnum       (1,1) {mustBePositive,mustBeInteger}
% [解説 L170] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                options.windowVector (:,1) double
% [解説 L171] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                options.windowName   char
            end
            
% [解説 L174] 条件「isfield(options,'windowVector')」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if isfield(options,'windowVector')
% [解説 L175] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.win」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                obj.win = options.windowVector;
% [解説 L176] この処理を実行して現在の計算状態を更新し、次の処理が必要とする状態へ進めます。
            elseif isfield(options,'windowName')
% [解説 L177] この関数の複数の計算結果を受け取り、obj、win、obj、diffWin にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                [obj.win,obj.diffWin] = DGTtool.getWindow(options.windowLength,options.windowName);
% [解説 L178] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L179] この関数の複数の計算結果を受け取り、obj、win、obj、diffWin にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                [obj.win,obj.diffWin] = DGTtool.getWindow(options.windowLength,'4termC5Nuttall');
            end
            
% [解説 L182] 条件「isfield(options,'FFTnum')」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if isfield(options,'FFTnum')
% [解説 L183] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.FFTnum」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
                obj.FFTnum = options.FFTnum;
% [解説 L184] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L185] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 obj.FFTnum」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
                obj.FFTnum = length(obj.win);
            end
            
% [解説 L188] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.shift」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            obj.shift  = options.windowShift;
% [解説 L189] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 obj.factorIdx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
            obj.factorIdx = struct('c',[],'d',[],'p',[],'q',[],'wIdx',[],'xIdx',[],'cIdx',[]);
        end
        
% [解説 L192] この関数は入力を受け取り、varargout を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function varargout = subsref(obj,data)
            %SUBSREF: Shortcut for DGT.
            %   Parenthesis notation (operator form) calls DGT.
            %
            %   Usage:
            %      X = F(x)
            %
            %   See also DGT
            
            % Note: F(x) is shortcut of subsref(F,struct('type','()','subs',{{x}}))
            
% [解説 L203] 値 data(1).type に応じて、入力形式・音源モデル・評価モードなどの処理経路を切り替えます。
            switch data(1).type
% [解説 L204] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                case '()' % Call DGT
                    % Many of validations are skipped for speed.
% [解説 L206] 条件「length(data) ~= 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                    if length(data) ~= 1
% [解説 L207] この条件では正しい計算ができないため、理由を示して停止します。
% [意図] 前提を満たさない入力で誤った数値結果を作らないよう、処理を停止するためです。
                        error 'DGT must be performed as F(x)'
                    end
% [解説 L209] 条件「numel(data.subs) ~= 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                    if numel(data.subs) ~= 1
% [解説 L210] 計算結果の解釈に影響する条件を利用者へ知らせます。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                        warning 'Only one input is allowed. Others will be ignored.'
                    end
                    
% [解説 L213] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                    x = data.subs{1};
                    
% [解説 L215] この関数の複数の計算結果を受け取り、varargout、nargout にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                    [varargout{1:nargout}] = DGT(obj,x);
                    
% [解説 L217] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                case '.' % implemented for access to properties and methods
% [解説 L218] 値 length(data) に応じて、入力形式・音源モデル・評価モードなどの処理経路を切り替えます。
                    switch length(data)
% [解説 L219] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                        case 1
% [解説 L220] この関数の複数の計算結果を受け取り、varargout、nargout にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                            [varargout{1:nargout}] = obj.(data.subs);
% [解説 L221] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                        case 2
% [解説 L222] この関数の複数の計算結果を受け取り、varargout、nargout にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                            [varargout{1:nargout}] = obj.(data(1).subs)(data(2).subs{:});
% [解説 L223] 直前までの条件に該当しなかった場合の代替処理へ進みます。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                        otherwise % not implemented for now
% [解説 L224] この条件では正しい計算ができないため、理由を示して停止します。
% [意図] 前提を満たさない入力で誤った数値結果を作らないよう、処理を停止するためです。
                            error 'Unexpected usage!'
                    end
                    
% [解説 L227] 直前までの条件に該当しなかった場合の代替処理へ進みます。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                otherwise % not implemented for now
% [解説 L228] この条件では正しい計算ができないため、理由を示して停止します。
% [意図] 前提を満たさない入力で誤った数値結果を作らないよう、処理を停止するためです。
                    error 'Unexpected usage!'
            end
        end
        
% [解説 L232] この関数は入力を受け取り、X,f,t を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function [X,f,t] = DGT(obj,x)
            %DGT: Compute spectrogram by DGT.
            %   Parenthesis notation can be used for shortcut.
            %   Normalized frequency f and sample index t can be returned.
            %
            %   Usage:
            %      X = F(x)
            %      [X,f,t] = F(x)
            %
            %   Note: The following notation gives the same result.
            %      X = F.DGT(x)
            %      X = DGT(F,x)
            %
            %   See also DGTtool, DGTtool/DGTtool, reassign
            
% [解説 L247] 条件「obj.FFTnum < length(obj.win)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if obj.FFTnum < length(obj.win)
% [解説 L248] 信号またはベクトルの長さを取得します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                x = zeroPadForFactorAlg(obj,x,length(x));
            end
            
% [解説 L251] 条件「size(x,1) < length(obj.win)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if size(x,1) < length(obj.win)
% [解説 L252] setFlag_winCalcInDGT を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                setFlag_winCalcInDGT(obj)
% [解説 L253] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 obj.win」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                obj.win = zeroPadForFactorAlg(obj,obj.win,length(obj.win));
% [解説 L254] 信号またはベクトルの長さを取得します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                x = zeroPadForFactorAlg(obj,x,length(obj.win));
            end
            
% [解説 L257] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 useFactorizationAlgorithm」へ代入します。
% [意図] 時間信号と時間周波数表現の軸・長さ・係数を整え、DGTと逆変換を対応させるためです。
            useFactorizationAlgorithm = size(x,1) == length(obj.win);
            
% [解説 L259] 条件「~isequal(obj.sigLen,size(x,1))」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isequal(obj.sigLen,size(x,1))
% [解説 L260] 条件「useFactorizationAlgorithm」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if useFactorizationAlgorithm
% [解説 L261] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 c」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                    c = lcm(obj.shift,obj.FFTnum);
% [解説 L262] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
                else
% [解説 L263] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 c」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                    c = obj.shift;
                end
% [解説 L265] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 obj.sigLen」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
                obj.sigLen = ceil(size(x,1)/c)*c;
            end
            
% [解説 L268] 条件「useFactorizationAlgorithm」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if useFactorizationAlgorithm
% [解説 L269] 条件「~isequal(obj.sigLen,size(x,1))」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if ~isequal(obj.sigLen,size(x,1))
% [解説 L270] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                    x = zeroPadForFactorAlg(obj,x,obj.sigLen);
                end
% [解説 L272] 条件「~isequal(obj.sigLen,length(obj.win))」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if ~isequal(obj.sigLen,length(obj.win))
% [解説 L273] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.win」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                    obj.win = zeroPadForFactorAlg(obj,obj.win,obj.sigLen);
                end
% [解説 L275] 条件「factorIdxMismatch(obj.shift,obj.FFTnum,size(x,1)/obj.shift,obj.factorIdx)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if factorIdxMismatch(obj.shift,obj.FFTnum,size(x,1)/obj.shift,obj.factorIdx)
% [解説 L276] computeIndexForFactorAlg を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます
                    computeIndexForFactorAlg(obj,size(x,1)/obj.shift,length(x),size(x,2));
                end
                
% [解説 L279] 条件「nargout > 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if nargout > 1
% [解説 L280] この関数の複数の計算結果を受け取り、X、f、t にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                    [X,f,t] = DGT_factorAlg(x,obj.win,obj.shift,obj.FFTnum,obj.factorIdx);
% [解説 L281] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
                else
% [解説 L282] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                    X = DGT_factorAlg(x,obj.win,obj.shift,obj.FFTnum,obj.factorIdx);
                end
% [解説 L284] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L285] 条件「nargout > 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if nargout > 1
% [解説 L286] この関数の複数の計算結果を受け取り、X、f、t にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                    [X,f,t] = DGT_usualAlg(x,obj.win,obj.shift,obj.FFTnum,obj.sigLen);
% [解説 L287] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
                else
% [解説 L288] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                    X = DGT_usualAlg(x,obj.win,obj.shift,obj.FFTnum,obj.sigLen);
                end
            end
        end
        
% [解説 L293] この関数は入力を受け取り、x を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function x = H(obj,X,sWin)
            %H: Inverse DGT.
            %   Synthesis window can be specified. This function is complex
            %   conjugate transpose of DGT if a window is not specified.
            %
            %   Usage:
            %      x = F.H(X)
            %      x = F.H(X,synthesisWindow)
            %
            %   See also pinv, makeWindowTight
            
            % Validation of input arguments is skipped for speed.
            % By default, this function uses analysis window.
% [解説 L306] 条件「~exist('sWin','var')」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~exist('sWin','var')
% [解説 L307] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 sWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                sWin = obj.win;
            end
            
% [解説 L310] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 signalLength」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
            signalLength = size(X,2)*obj.shift;
% [解説 L311] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 useFactorizationAlgorithm」へ代入します。
% [意図] 時間信号と時間周波数表現の軸・長さ・係数を整え、DGTと逆変換を対応させるためです。
            useFactorizationAlgorithm = signalLength <= length(sWin);
            
% [解説 L313] 条件「useFactorizationAlgorithm」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if useFactorizationAlgorithm
% [解説 L314] 条件「length(sWin) ~= length(obj.win)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if length(sWin) ~= length(obj.win)
% [解説 L315] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                    X = changeDGTdef(obj,X);
                end
% [解説 L317] 条件「length(sWin) ~= signalLength」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if length(sWin) ~= signalLength
% [解説 L318] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 sWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                    sWin = zeroPadForFactorAlg(obj,sWin,length(sWin));
% [解説 L319] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                    X = [X, X(:,1:(length(sWin)/obj.shift)-size(X,2),:)];
                end
% [解説 L321] 条件「isempty(obj.factorIdx) || sizeMismatch(obj.factorIdx,X) || factorIdxMismatch(obj.shift,obj.FFTnum,size(X,2),obj.factorIdx)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if isempty(obj.factorIdx) || sizeMismatch(obj.factorIdx,X) || factorIdxMismatch(obj.shift,obj.FFTnum,size(X,2),obj.factorIdx)
% [解説 L322] computeIndexForFactorAlg を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます
                    computeIndexForFactorAlg(obj,size(X,2),length(sWin),size(X,3));
                end
                
% [解説 L325] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                x = invDGT_factorAlg(X,sWin,obj.shift,obj.FFTnum,obj.factorIdx);
% [解説 L326] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L327] 条件「~isequal(size(obj.OLAindex,1:3),[length(sWin) size(X,[2 3])])」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if ~isequal(size(obj.OLAindex,1:3),[length(sWin) size(X,[2 3])])
% [解説 L328] computeIndexForOLA を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                    computeIndexForOLA(obj,X,sWin);
                end
                
% [解説 L331] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                x = invDGT_usualAlg(X,sWin,obj.FFTnum,obj.OLAindex);
            end
            
% [解説 L334] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 obj.sigLen」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
            obj.sigLen = size(x,1);
        end
        
% [解説 L337] この関数は入力を受け取り、x を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function x = pinv(obj,X)
            %PINV: Inverse DGT with perfect reconstruction.
            %   Signal is reconstructed using canonical dual window.
            %
            %   Usage:
            %      x = F.pinv(X)
            %
            %   See also DGTtool, DGT, H, computeCanonicalDual
            
% [解説 L346] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 signalLength」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
            signalLength = size(X,2)*obj.shift;
            
% [解説 L348] 条件「isempty(obj.dualWin) || ~obj.isDual || signalLength < length(obj.dualWin)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if isempty(obj.dualWin) || ~obj.isDual || signalLength < length(obj.dualWin)
% [解説 L349] 候補の最大値を取り、ピークや最良候補を求めます。結果を「処理途中で使う変数 newLen」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
                newLen = max(signalLength,length(obj.dualWin));
% [解説 L350] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 newLen」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
                newLen = DGTtool.sigLenForFactorDGT(newLen,obj.shift,obj.FFTnum);
% [解説 L351] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.dualWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                obj.dualWin = getCanonicalDualWin(obj,newLen);
            end
% [解説 L353] 条件「~obj.isCanonical」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~obj.isCanonical
% [解説 L354] 計算結果の解釈に影響する条件を利用者へ知らせます。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                warning 'This is not pseudoinverse because the given dual window is not canonical. Consider using F.H(x,F.dualWin) for non-canonical dual.'
            end
            
% [解説 L357] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 useFactorizationAlgorithm」へ代入します。
% [意図] 時間信号と時間周波数表現の軸・長さ・係数を整え、DGTと逆変換を対応させるためです。
            useFactorizationAlgorithm = signalLength <= length(obj.dualWin);
            
% [解説 L359] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
            x = H(obj,X,obj.dualWin);
            
% [解説 L361] 条件「useFactorizationAlgorithm && length(obj.win) ~= length(obj.dualWin)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if useFactorizationAlgorithm && length(obj.win) ~= length(obj.dualWin)
% [解説 L362] 信号またはベクトルの長さを取得します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                x = circshift(x,-length(obj.win)+obj.shift);
            end
        end
        
% [解説 L366] この関数は入力を受け取り、reassignedS,f,t,X,IF,GD,dXdt,dXdf を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function [reassignedS,f,t,X,IF,GD,dXdt,dXdf] = reassign(obj,x,epsilon)
            %REASSIGN: Compute reassigned spectrogram.
            %   Reassigned spectrogram is a sparse time-frequency representation.
            %
            %   Usage:
            %      rS = F.reassign(x)
            %      rS = F.reassign(x,epsilon)
            %      [rS,f,t,X,IF,GD,dXdt,dXdf] = F.reassign(x)
            %
            %   Input:
            %      x       - Time-domain signal (column vectors)
            %      epsilon - Small number used for avoiding zero-division
            %
            %   Output:
            %      rS   - Reassigned spectrogram (non-negative valued)
            %      f    - Normalized frequency (column vector)
            %      t    - Time indices (row vector)
            %      X    - Spectrogram (complex-valued)
            %      IF   - Instantaneous frequency (time-derivative of phase)
            %      GD   - Group delay (frequency-derivative of phase)
            %      dXdt - Time-derivative of spectrogram
            %      dXdf - Frequency-derivative of spectrogram
            %
            %   See also plotReassign
% [解説 L390] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L391] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                obj
% [解説 L392] 複数音源が混ざった観測時間信号 を (:,:)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
                x       (:,:) double
% [解説 L393] 処理途中で使う変数 epsilon を (1,1)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
                epsilon (1,1) double {mustBePositive} = 1e-12;
            end
            
% [解説 L396] この関数の複数の計算結果を受け取り、maxIdx にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
            [~,maxIdx] = max(obj.win);
% [解説 L397] 条件「(maxIdx - length(obj.win)/2) > 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if (maxIdx - length(obj.win)/2) > 1
% [解説 L398] 計算結果の解釈に影響する条件を利用者へ知らせます。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                warning 'Window seems improper. Consider using F.setWindow.'
            end
            
% [解説 L401] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 nWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
            nWin = obj.win / sum(obj.win);
            
% [解説 L403] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 cg」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            cg = centerOfGravity(nWin);
% [解説 L404] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 tRamp」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            tRamp = (-length(nWin)/2:length(nWin)/2-1)' + mod(length(nWin),2)/2; % mod for odd/even cases
% [解説 L405] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 tRamp」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            tRamp = circshift(tRamp,floor(cg)-floor(length(nWin)/2)-1);
% [解説 L406] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 tWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
            tWin = tRamp .* nWin;
            
% [解説 L408] 条件「isempty(obj.diffWin)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if isempty(obj.diffWin)
% [解説 L409] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 obj.sigLen」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
                obj.sigLen = DGTtool.sigLenForFactorDGT(size(x,1),obj.shift,obj.FFTnum);
% [解説 L410] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 dWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                dWin = DGTtool.computeNumericalDiffWin(nWin,obj.sigLen) * obj.FFTnum / length(obj.win);
                
% [解説 L412] 時間信号と時間周波数表現の間を変換します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                x    = DGTtool.extendSignalByZeroPad(x,obj.sigLen);
% [解説 L413] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 tWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                tWin = DGTtool.extendSignalByZeroPad(tWin,obj.sigLen);
% [解説 L414] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 nWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                nWin = DGTtool.extendSignalByZeroPad(nWin,obj.sigLen);
% [解説 L415] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L416] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 obj.sigLen」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
                obj.sigLen = ceil(size(x,1)/obj.shift) * obj.shift;
% [解説 L417] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 dWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                dWin = obj.diffWin * obj.FFTnum / length(obj.win) / sum(obj.win);
            end
            
% [解説 L420] 条件「size(x,1) < length(nWin) || size(x,1) < length(dWin)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if size(x,1) < length(nWin) || size(x,1) < length(dWin)
% [解説 L421] この条件では正しい計算ができないため、理由を示して停止します。
% [意図] 前提を満たさない入力で誤った数値結果を作らないよう、処理を停止するためです。
                error 'Signal is shorter than window.'
            end
% [解説 L423] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 useFactorizationAlgorithm」へ代入します。
% [意図] 時間信号と時間周波数表現の軸・長さ・係数を整え、DGTと逆変換を対応させるためです。
            useFactorizationAlgorithm = obj.sigLen == length(dWin);
            
% [解説 L425] 条件「useFactorizationAlgorithm」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if useFactorizationAlgorithm
% [解説 L426] 条件「factorIdxMismatch(obj.shift,obj.FFTnum,size(x,1)/obj.shift,obj.factorIdx)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if factorIdxMismatch(obj.shift,obj.FFTnum,size(x,1)/obj.shift,obj.factorIdx)
% [解説 L427] computeIndexForFactorAlg を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます
                    computeIndexForFactorAlg(obj,size(x,1)/obj.shift,length(x),size(x,2));
                end
% [解説 L429] この関数の複数の計算結果を受け取り、X、f、t にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                [X,f,t] = DGT_factorAlg(x,nWin,obj.shift,obj.FFTnum,obj.factorIdx);
% [解説 L430] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 dXdt」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                dXdt = DGT_factorAlg(x,dWin,obj.shift,obj.FFTnum,obj.factorIdx);
% [解説 L431] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 dXdf」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                dXdf = DGT_factorAlg(x,tWin,obj.shift,obj.FFTnum,obj.factorIdx);
% [解説 L432] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L433] この関数の複数の計算結果を受け取り、X、f、t にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                [X,f,t] = DGT_usualAlg(x,nWin,obj.shift,obj.FFTnum,obj.sigLen);
% [解説 L434] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 dXdt」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                dXdt = DGT_usualAlg(x,dWin,obj.shift,obj.FFTnum,obj.sigLen);
% [解説 L435] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 dXdf」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                dXdf = DGT_usualAlg(x,tWin,obj.shift,obj.FFTnum,obj.sigLen);
            end
            
% [解説 L438] 複素値または符号付き値を振幅・大きさへ変換します。結果を「正解音源または正解音像の時間信号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            S = abs(X).^2;
% [解説 L439] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 Splus」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            Splus = S + epsilon;
            
% [解説 L441] 条件「= -imag(dXdt.*conj(X)./Splus);」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            IF = -imag(dXdt.*conj(X)./Splus);
% [解説 L442] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 GD」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            GD =  real(dXdf.*conj(X)./Splus);
            
% [解説 L444] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 fNew」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            fNew = f*obj.FFTnum + IF;
% [解説 L445] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 tNew」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            tNew = t + GD;
            
% [解説 L447] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 fNewIdx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
            fNewIdx = round(fNew(:)) + 1;
% [解説 L448] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 tNewIdx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
            tNewIdx = round(tNew(:)/obj.shift);
% [解説 L449] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 tNewIdx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
            tNewIdx = mod(tNewIdx,size(X,2)) + 1;
            
% [解説 L451] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 idx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
            idx = (fNewIdx >= 1) & (fNewIdx <= size(X,1));
% [解説 L452] 要素数を保ったままshapeを組み替えます。結果を「処理途中で使う変数 chIdx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
            chIdx = ones(size(X)) .* reshape(1:size(X,3),1,1,[]);
% [解説 L453] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 reassignedS」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
            reassignedS = accumarray([fNewIdx(idx) tNewIdx(idx) chIdx(idx)], S(idx), size(X,1:3));
            
% [解説 L455] 条件「nargin > 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if nargin > 1
% [解説 L456] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                t = mod(t,size(X,2)*obj.shift);
% [解説 L457] この関数の複数の計算結果を受け取り、tRot にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                [~,tRot] = min(t);
% [解説 L458] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                X = circshift(X,-tRot+1,2);
% [解説 L459] 条件「= circshift(IF,-tRot+1,2);」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                IF = circshift(IF,-tRot+1,2);
% [解説 L460] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 GD」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                GD = circshift(GD,-tRot+1,2);
% [解説 L461] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 dXdt」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                dXdt = circshift(dXdt,-tRot+1,2);
% [解説 L462] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 dXdf」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                dXdf = circshift(dXdf,-tRot+1,2);
% [解説 L463] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                t = circshift(t,-tRot+1);
            end
        end
        
% [解説 L467] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function setWindow(obj,windowLength,windowName,varargin)
            %setWindow: Set a new window (current window is deleted).
            %   All input variables are directly passed to DGTtool.getWindow.
            %
            %   Usage:
            %      F.setWindow(windowLength,windowName)
            %      F.setWindow(___,Name,Value)
            %
            %   See also getWindow, plotWin
            
% [解説 L477] この関数の複数の計算結果を受け取り、obj、win、obj、diffWin にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
            [obj.win,obj.diffWin] = DGTtool.getWindow(windowLength,windowName,varargin{:});
        end
        
% [解説 L480] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function makeWindowTight(obj,signalLength)
            %makeWindowTight: Replace window by its canonical tight window.
            %   Canonical tight window makes F.H(X) and F.pinv(X) equal.
            %   Signal length must be specified whenever FFTnum < winLen.
            %
            %   Usage:
            %      F.makeWindowTight
            %      F.makeWindowTight(signalLength)
            %
            %   See also computeCanonicalTight, plotWin
            
% [解説 L491] 条件「~exist('signalLength','var')」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~exist('signalLength','var')
% [解説 L492] 候補の最大値を取り、ピークや最良候補を求めます。結果を「処理途中で使う変数 signalLength」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
                signalLength = max(obj.sigLen,length(obj.win));
% [解説 L493] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L494] 候補の最大値を取り、ピークや最良候補を求めます。結果を「処理途中で使う変数 signalLength」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
                signalLength = max(signalLength,length(obj.win));
% [解説 L495] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.sigLen」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
                obj.sigLen = ceil(signalLength/obj.shift)*obj.shift;
            end
% [解説 L497] 条件「isempty(signalLength)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if isempty(signalLength)
% [解説 L498] 条件「obj.FFTnum < length(obj.win)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if obj.FFTnum < length(obj.win)
% [解説 L499] この条件では正しい計算ができないため、理由を示して停止します。
% [意図] 前提を満たさない入力で誤った数値結果を作らないよう、処理を停止するためです。
                    error 'Must specify signal length, e.g., F.makeWindowTight(length(x)).'
% [解説 L500] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
                else
% [解説 L501] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 signalLength」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
                    signalLength = length(obj.win);
                end
            end
% [解説 L504] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 obj.win」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
            obj.win = DGTtool.computeCanonicalTight(obj.win,obj.shift,obj.FFTnum,signalLength);
% [解説 L505] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.dualWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
            obj.dualWin = obj.win;
        end
        
% [解説 L508] この関数は入力を受け取り、X を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function X = changeDGTdef(obj,X)
            %changeDGTdef: Phase is modified to change definition of DGT.
            %
            %   Usage:
            %      X = F.changeDGTdef(X)
            %
            %   See also undoChangeDGTdef
            
% [解説 L516] 条件「~isequal(size(obj.defConverter),size(X,1:2))」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isequal(size(obj.defConverter),size(X,1:2))
% [解説 L517] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.defConverter」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                obj.defConverter = calculateDefConverter(X,obj.shift,obj.FFTnum);
            end
% [解説 L519] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
            X = X .* obj.defConverter;
        end
        
% [解説 L522] この関数は入力を受け取り、X を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function X = undoChangeDGTdef(obj,X)
            %undoChangeDGTdef: Reset phase modified by changeDGTdef.
            %
            %   Usage:
            %      X = F.undoChangeDGTdef(X)
            %
            %   See also changeDGTdef
            
% [解説 L530] 条件「~isequal(size(obj.defConverter),size(X,1:2))」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isequal(size(obj.defConverter),size(X,1:2))
% [解説 L531] この条件では正しい計算ができないため、理由を示して停止します。
% [意図] 前提を満たさない入力で誤った数値結果を作らないよう、処理を停止するためです。
                error 'Parameters seem different.'
            end
% [解説 L533] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
            X = X .* conj(obj.defConverter);
        end
        
% [解説 L536] この関数は入力を受け取り、X を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function X = makeZeroPhase(obj,X)
            %makeZeroPhase: Phase is modified to remove linear phase component of window.
            %
            %   Usage:
            %      X = F.makeZeroPhase(X)
            %
            %   See also undoMakeZeroPhase
            
% [解説 L544] 条件「~isequal(size(obj.zeroPhaseConverter,1),size(X,1))」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isequal(size(obj.zeroPhaseConverter,1),size(X,1))
% [解説 L545] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.zeroPhaseConverter」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                obj.zeroPhaseConverter = calculateZeroPhaseConverter(X,obj.FFTnum,obj.win);
            end
% [解説 L547] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
            X = X .* obj.zeroPhaseConverter;
        end
        
% [解説 L550] この関数は入力を受け取り、X を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function X = undoMakeZeroPhase(obj,X)
            %undoMakeZeroPhase: Reset phase modified by makeZeroPhase.
            %
            %   Usage:
            %      X = F.undoMakeZeroPhase(X)
            %
            %   See also makeZeroPhase
            
% [解説 L558] 条件「~isequal(size(obj.zeroPhaseConverter,1),size(X,1))」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isequal(size(obj.zeroPhaseConverter,1),size(X,1))
% [解説 L559] この条件では正しい計算ができないため、理由を示して停止します。
% [意図] 前提を満たさない入力で誤った数値結果を作らないよう、処理を停止するためです。
                error 'Parameters seem different.'
            end
% [解説 L561] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
            X = X .* conj(obj.zeroPhaseConverter);
        end
    end
    
% [解説 L565] DGT、逆変換、描画、窓生成、設定更新を行うメソッド群を定義します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
    methods % plot
        
% [解説 L567] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function plot(obj,x,fs,options)
            %PLOT: Compute spectrogram and display it.
            %   Signal and spectrum are also displayed.
            %
            %   Usage:
            %      F.plot(x)
            %      F.plot(x,fs)
            %      F.plot(___,Name,Value)
            %
            %   Options:
            %      'range'     (number, default = 80 [dB])
            %      'trunc'     (number, default = 0  [dB])
            %      'normalize' (true/false, default = true)
            %
            %   See also plotPhase, plotReassign, DGT
% [解説 L582] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L583] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                obj
% [解説 L584] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                x   (:,:) {mustBeNumeric,mustBeSkinny}
% [解説 L585] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                fs  (1,1) {mustBePositive} = 1
% [解説 L586] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 options.range」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                options.range = 80
% [解説 L587] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 options.trunc」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                options.trunc = 0
% [解説 L588] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 options.normalize」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                options.normalize = true
            end
            
% [解説 L591] 条件「options.normalize」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if options.normalize
% [解説 L592] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 normConst」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                normConst = sum(obj.win) / 2;
% [解説 L593] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L594] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 normConst」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                normConst = 1;
            end
            
% [解説 L597] この関数の複数の計算結果を受け取り、X、f、t にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
            [X,f,t] = DGT(obj,x);
% [解説 L598] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
            X = X / normConst;
            
% [解説 L600] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            t = mod(t,size(X,2)*obj.shift);
% [解説 L601] この関数の複数の計算結果を受け取り、tRot にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
            [~,tRot] = min(t);
% [解説 L602] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
            X = circshift(X,-tRot+1,2);
% [解説 L603] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            t = circshift(t,-tRot+1);
            
% [解説 L605] 条件「fs == 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if fs == 1
% [解説 L606] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 fsPlot」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                fsPlot = fs;
% [解説 L607] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 fUnit」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                fUnit = '[periods/sample]';
% [解説 L608] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 tUnit」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                tUnit = '[samples]';
% [解説 L609] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L610] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 fsPlot」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                fsPlot = fs/1000;
% [解説 L611] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「現在処理している周波数ビン番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                f = f*fsPlot;
% [解説 L612] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                t = t/fs;
% [解説 L613] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 fUnit」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                fUnit = '[kHz]';
% [解説 L614] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 tUnit」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                tUnit = '[s]';
            end
            
% [解説 L617] 処理途中で使う変数 n を 1:size(X,3) の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
            for n = 1:size(X,3)
% [解説 L618] 複素値または符号付き値を振幅・大きさへ変換します。結果を「正解音源または正解音像の時間信号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                s = 20*log10(abs(fft(x(:,n)/normConst)));
% [解説 L619] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 xLonger」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                xLonger = buffer(x(:,n),size(X,2)*obj.shift);
                
% [解説 L621] 信号・スペクトログラム・コスト推移を図にし、分離状態や収束挙動を目で確認できるようにします。
                figure
% [解説 L622] tiledlayout を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します
                tiledlayout(10,14,'TileSpacing','none','Padding','compact')
                
% [解説 L624] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 ax_s」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                ax_s = nexttile(1,[8 2]);
% [解説 L625] 信号・スペクトログラム・コスト推移を図にし、分離状態や収束挙動を目で確認できるようにします。
                plot(s(1:floor(length(s)/2)+1),(0:floor(length(s)/2))/length(s)*fsPlot)
% [解説 L626] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                xlim(max(s)-[options.range 0]-options.trunc)
% [解説 L627] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                ylim([0 floor(length(s)/2)/length(s)*fsPlot])
% [解説 L628] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                set(gca,'FontSize',10,'xDir','reverse')
% [解説 L629] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                ylabel(['Frequency ' fUnit],'FontSize',12)
                
% [解説 L631] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 ax_X」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                ax_X = nexttile(3,[8 12]);
% [解説 L632] 信号・スペクトログラム・コスト推移を図にし、分離状態や収束挙動を目で確認できるようにします。
                imagesc(t,f,20*log10(abs(X(:,:,n))))
% [解説 L633] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                axis xy off
% [解説 L634] caxis を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                caxis(max(caxis)-[options.range 0]-options.trunc)
% [解説 L635] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                xlim([0 size(x,1)-1]/fs)
% [解説 L636] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                ylim([0 floor(length(s)/2)/length(s)*fsPlot])
                
% [解説 L638] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 ax_x」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                ax_x = nexttile(115,[2 12]);
% [解説 L639] 信号・スペクトログラム・コスト推移を図にし、分離状態や収束挙動を目で確認できるようにします。
                plot((0:length(xLonger)-1)/fs,xLonger)
% [解説 L640] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                axis tight
% [解説 L641] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                xlim([0 size(x,1)-1]/fs)
% [解説 L642] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                set(gca,'FontSize',10)
% [解説 L643] yticks を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                yticks(0)
% [解説 L644] yline を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                yline(0)
% [解説 L645] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                xlabel(['Time ' tUnit],'FontSize',12)
                
% [解説 L647] nexttile を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                nexttile(127,[1 1])
% [解説 L648] image を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 要素数を保ったままshapeを組み替えます
                image([-1 1],ax_X.CLim,reshape(colormap,[],1,3))
% [解説 L649] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                axis xy
% [解説 L650] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                set(gca,'FontSize',10)
% [解説 L651] xticks を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                xticks([])
% [解説 L652] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                ylabel('Power [dB]','FontSize',11)
                
% [解説 L654] linkaxes を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します
                linkaxes([ax_X ax_x],'x')
% [解説 L655] linkaxes を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します
                linkaxes([ax_X ax_s],'y')
            end
        end
        
% [解説 L659] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function plotPhase(obj,x,fs,options)
            %plotPhase: Visualize phase.
            %   Color and brightness represent phase and magnitude, respectively.
            %   F.makeZeroPhase is applied for better visibility.
            %
            %   Usage:
            %      F.plotPhase(x)
            %      F.plotPhase(x,fs)
            %      F.plotPhase(___,Name,Value)
            %
            %   Options:
            %      'range'     (number, default = 40 [dB])
            %      'trunc'     (number, default = 15 [dB])
            %      'normalize' (true/false, default = true)
            %
            %   See also plot, plotReassign, makeZeroPhase
% [解説 L675] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L676] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                obj
% [解説 L677] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                x   (:,:) {mustBeNumeric,mustBeSkinny}
% [解説 L678] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                fs  (1,1) {mustBePositive} = 1
% [解説 L679] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 options.range」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                options.range = 40
% [解説 L680] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 options.trunc」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                options.trunc = 15
% [解説 L681] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 options.normalize」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                options.normalize = true
            end
            
% [解説 L684] 条件「options.normalize」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if options.normalize
% [解説 L685] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 normConst」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                normConst = sum(obj.win) / 2;
% [解説 L686] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L687] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 normConst」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                normConst = 1;
            end
            
% [解説 L690] この関数の複数の計算結果を受け取り、X、f、t にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
            [X,f,t] = DGT(obj,x);
% [解説 L691] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
            X = X / normConst;
% [解説 L692] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
            X = makeZeroPhase(obj,X);
            
% [解説 L694] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            t = mod(t,size(X,2)*obj.shift);
% [解説 L695] この関数の複数の計算結果を受け取り、tRot にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
            [~,tRot] = min(t);
% [解説 L696] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
            X = circshift(X,-tRot+1,2);
% [解説 L697] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            t = circshift(t,-tRot+1);
            
% [解説 L699] 条件「fs == 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if fs == 1
% [解説 L700] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 fsPlot」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                fsPlot = fs;
% [解説 L701] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 fUnit」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                fUnit = '[periods/sample]';
% [解説 L702] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 tUnit」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                tUnit = '[samples]';
% [解説 L703] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L704] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 fsPlot」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                fsPlot = fs/1000;
% [解説 L705] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「現在処理している周波数ビン番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                f = f*fsPlot;
% [解説 L706] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                t = t/fs;
% [解説 L707] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 fUnit」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                fUnit = '[kHz]';
% [解説 L708] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 tUnit」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                tUnit = '[s]';
            end
            
% [解説 L711] 処理途中で使う変数 n を 1:size(X,3) の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
            for n = 1:size(X,3)
% [解説 L712] 複素値または符号付き値を振幅・大きさへ変換します。結果を「処理途中で使う変数 A」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                A = 20*log10(abs(X(:,:,n)));
% [解説 L713] 候補の最大値を取り、ピークや最良候補を求めます。結果を「処理途中で使う変数 maxC」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                maxC = max(A(:)) - options.trunc;
% [解説 L714] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 minC」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                minC = maxC - options.range;
                
% [解説 L716] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 A」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                A = rescale(A,'InputMin',minC,'InputMax',maxC);
% [解説 L717] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 P」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                P = rescale(angle(X(:,:,n)),'InputMin',-pi,'InputMax',pi);
% [解説 L718] 初期値または重みとして使う1の配列を作ります。結果を「処理途中で使う変数 C」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                C = hsv2rgb(cat(3,P,ones(size(X(:,:,n))),A));
                
% [解説 L720] 複素値または符号付き値を振幅・大きさへ変換します。結果を「正解音源または正解音像の時間信号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                s = 20*log10(abs(fft(x(:,n)/normConst)));
% [解説 L721] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 xLonger」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                xLonger = buffer(x(:,n),size(X,2)*obj.shift);
                
% [解説 L723] 信号・スペクトログラム・コスト推移を図にし、分離状態や収束挙動を目で確認できるようにします。
                figure
% [解説 L724] tiledlayout を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します
                tiledlayout(10,14,'TileSpacing','none','Padding','compact')
                
% [解説 L726] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 ax_s」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                ax_s = nexttile(1,[8 2]);
% [解説 L727] 信号・スペクトログラム・コスト推移を図にし、分離状態や収束挙動を目で確認できるようにします。
                plot(s(1:floor(length(s)/2)+1),(0:floor(length(s)/2))/length(s)*fsPlot)
% [解説 L728] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                xlim(max(s)-[options.range 0]-options.trunc)
% [解説 L729] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                ylim([0 floor(length(s)/2)/length(s)*fsPlot])
% [解説 L730] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                set(gca,'FontSize',10,'xDir','reverse')
% [解説 L731] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                ylabel(['Frequency ' fUnit],'FontSize',12)
                
% [解説 L733] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 ax_X」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                ax_X = nexttile(3,[8 12]);
% [解説 L734] image を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                image(t,f,C)
% [解説 L735] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                axis xy off
% [解説 L736] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                xlim([0 size(x,1)-1]/fs)
% [解説 L737] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                ylim([0 floor(length(s)/2)/length(s)*fsPlot])
                
% [解説 L739] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 ax_x」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                ax_x = nexttile(115,[2 12]);
% [解説 L740] 信号・スペクトログラム・コスト推移を図にし、分離状態や収束挙動を目で確認できるようにします。
                plot((0:length(xLonger)-1)/fs,xLonger)
% [解説 L741] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                axis tight
% [解説 L742] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                xlim([0 size(x,1)-1]/fs)
% [解説 L743] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                set(gca,'FontSize',10)
% [解説 L744] yticks を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                yticks(0)
% [解説 L745] yline を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                yline(0)
% [解説 L746] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                xlabel(['Time ' tUnit],'FontSize',12)
                
% [解説 L748] nexttile を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                nexttile(127,[1 1])
% [解説 L749] image を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 初期値または重みとして使う1の配列を作ります
                image([-1 1],[minC maxC],hsv2rgb(cat(3,repmat(linspace(0,1,64),128,1),ones(128,64),repmat(linspace(0,1,128)',1,64))))
% [解説 L750] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                axis xy
% [解説 L751] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                set(gca,'FontSize',10)
% [解説 L752] xticks を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                xticks(-1:1)
% [解説 L753] xticklabels を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 逆行列を明示せず線形方程式を解き、数値的に安定な係数を求めます
                xticklabels({'-\pi',0,'\pi'})
% [解説 L754] xtickangle を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                xtickangle(0)
% [解説 L755] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                xlabel('Phase [rad]','FontSize',11)
% [解説 L756] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                ylabel('Power [dB]','FontSize',11)
                
% [解説 L758] linkaxes を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します
                linkaxes([ax_X ax_x],'x')
% [解説 L759] linkaxes を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します
                linkaxes([ax_X ax_s],'y')
            end
        end
        
% [解説 L763] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function plotReassign(obj,x,fs,options)
            %plotReassign: Compute reassigned spectrogram and display it.
            %   Resolution of reassigned spectrogram depends on shift and FFTnum.
            %
            %   Usage:
            %      F.plotReassign(x)
            %      F.plotReassign(x,fs)
            %      F.plotReassign(___,Name,Value)
            %
            %   Options:
            %      'range'     (number, default = 100 [dB])
            %      'trunc'     (number, default = 0   [dB])
            %      'normalize' (true/false, default = true)
            %
            %   See also plot, plotPhase, reassign
% [解説 L778] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L779] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                obj
% [解説 L780] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                x   (:,:) {mustBeNumeric,mustBeSkinny}
% [解説 L781] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                fs  (1,1) {mustBePositive} = 1
% [解説 L782] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 options.range」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                options.range = 100
% [解説 L783] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 options.trunc」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                options.trunc = 0
% [解説 L784] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 options.normalize」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                options.normalize = true
            end
            
% [解説 L787] この関数の複数の計算結果を受け取り、X、f、t にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
            [X,f,t] = reassign(obj,x);
            
% [解説 L789] 条件「options.normalize」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if options.normalize
% [解説 L790] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 normConst」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                normConst = sum(obj.win);
% [解説 L791] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L792] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 normConst」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                normConst = 1;
% [解説 L793] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                X = X * sum(obj.win);
            end
            
% [解説 L796] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            t = mod(t,size(X,2)*obj.shift);
% [解説 L797] この関数の複数の計算結果を受け取り、tRot にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
            [~,tRot] = min(t);
% [解説 L798] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
            X = circshift(X,-tRot+1,2);
% [解説 L799] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            t = circshift(t,-tRot+1);
            
% [解説 L801] 条件「fs == 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if fs == 1
% [解説 L802] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 fsPlot」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                fsPlot = fs;
% [解説 L803] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 fUnit」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                fUnit = '[periods/sample]';
% [解説 L804] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 tUnit」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                tUnit = '[samples]';
% [解説 L805] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L806] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 fsPlot」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                fsPlot = fs/1000;
% [解説 L807] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「現在処理している周波数ビン番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                f = f*fsPlot;
% [解説 L808] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                t = t/fs;
% [解説 L809] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 fUnit」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                fUnit = '[kHz]';
% [解説 L810] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 tUnit」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                tUnit = '[s]';
            end
            
% [解説 L813] 処理途中で使う変数 n を 1:size(X,3) の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
            for n = 1:size(X,3)
% [解説 L814] 複素値または符号付き値を振幅・大きさへ変換します。結果を「正解音源または正解音像の時間信号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                s = 20*log10(abs(fft(x(:,n)/normConst)));
% [解説 L815] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 xLonger」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
                xLonger = buffer(x(:,n),size(X,2)*obj.shift);
                
% [解説 L817] 信号・スペクトログラム・コスト推移を図にし、分離状態や収束挙動を目で確認できるようにします。
                figure
% [解説 L818] tiledlayout を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します
                tiledlayout(10,14,'TileSpacing','none','Padding','compact')
                
% [解説 L820] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 ax_s」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                ax_s = nexttile(1,[8 2]);
% [解説 L821] 信号・スペクトログラム・コスト推移を図にし、分離状態や収束挙動を目で確認できるようにします。
                plot(s(1:floor(length(s)/2)+1),(0:floor(length(s)/2))/length(s)*fsPlot)
% [解説 L822] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                xlim(max(s)-[options.range 0]-options.trunc)
% [解説 L823] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                ylim([0 floor(length(s)/2)/length(s)*fsPlot])
% [解説 L824] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                set(gca,'FontSize',10,'xDir','reverse')
% [解説 L825] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                ylabel(['Frequency ' fUnit],'FontSize',12)
                
% [解説 L827] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 ax_X」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                ax_X = nexttile(3,[8 12]);
% [解説 L828] 信号・スペクトログラム・コスト推移を図にし、分離状態や収束挙動を目で確認できるようにします。
                imagesc(t,f,20*log10(abs(X(:,:,n))))
% [解説 L829] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                axis xy off
% [解説 L830] caxis を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                caxis(max(caxis)-[options.range 0]-options.trunc)
% [解説 L831] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                xlim([0 size(x,1)-1]/fs)
% [解説 L832] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                ylim([0 floor(length(s)/2)/length(s)*fsPlot])
                
% [解説 L834] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 ax_x」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                ax_x = nexttile(115,[2 12]);
% [解説 L835] 信号・スペクトログラム・コスト推移を図にし、分離状態や収束挙動を目で確認できるようにします。
                plot((0:length(xLonger)-1)/fs,xLonger)
% [解説 L836] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                axis tight
% [解説 L837] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                xlim([0 size(x,1)-1]/fs)
% [解説 L838] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                set(gca,'FontSize',10)
% [解説 L839] yticks を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                yticks(0)
% [解説 L840] yline を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                yline(0)
% [解説 L841] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                xlabel(['Time ' tUnit],'FontSize',12)
                
% [解説 L843] nexttile を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                nexttile(127,[1 1])
% [解説 L844] image を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 要素数を保ったままshapeを組み替えます
                image([-1 1],ax_X.CLim,reshape(colormap,[],1,3))
% [解説 L845] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                axis xy
% [解説 L846] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                set(gca,'FontSize',10)
% [解説 L847] xticks を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                xticks([])
% [解説 L848] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
                ylabel('Power [dB]','FontSize',11)
                
% [解説 L850] linkaxes を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します
                linkaxes([ax_X ax_x],'x')
% [解説 L851] linkaxes を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します
                linkaxes([ax_X ax_s],'y')
            end
        end
        
% [解説 L855] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function plotWin(obj)
            %plotWin: Display currently available windows.
            %   Some windows may not appear if they are not calculated yet.
            %
            %   Usage:
            %      F.plotWin
            %
            %   See also setWindow, getWindow
            
% [解説 L864] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 h」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            h = figure;
% [解説 L865] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 h.Position(4)」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            h.Position(4) = h.Position(4)/2;
% [解説 L866] tiledlayout を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します
            tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
            
% [解説 L868] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
            nexttile
% [解説 L869] windowStylePlot を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
            windowStylePlot(obj.win)
% [解説 L870] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
            title('window','FontSize',12)
            
% [解説 L872] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
            nexttile
% [解説 L873] windowStylePlot を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
            windowStylePlot(obj.dualWin)
% [解説 L874] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
            title('dual window','FontSize',12)
            
% [解説 L876] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
            nexttile
% [解説 L877] windowStylePlot を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
            windowStylePlot(obj.diffWin)
% [解説 L878] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
            title('differential window','FontSize',12)
        end
        
    end
    
% [解説 L883] DGT、逆変換、描画、窓生成、設定更新を行うメソッド群を定義します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
    methods (Static)
        
% [解説 L885] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        function winList = windowList
            %winList: List of window names acceptable in DGTtool.
            %
            %   Usage:
            %      c = DGTtool.windowList
            %
            %   See also getWindow, DGTtool, DGTtool/DGTtool
            
% [解説 L893] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 winList」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
            winList = {
% [解説 L894] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                'ForKitamuraHamming'
% [解説 L895] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                'Hann'
% [解説 L896] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                'Blackman'
% [解説 L897] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                '3termC1Nuttall'
% [解説 L898] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                '3termC3Nuttall'
% [解説 L899] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                '4termC1Nuttall'
% [解説 L900] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                '4termC3Nuttall'
% [解説 L901] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                '4termC5Nuttall'
% [解説 L902] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                'Gauss'
% [解説 L903] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                'Slepian'
% [解説 L904] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                'Chebyshev'
% [解説 L905] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                };
        end
        
% [解説 L908] この関数は入力を受け取り、win,diffWin を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function [win,diffWin] = getWindow(windowLength,windowName,options)
            %getWindow: Compute specified window.
            %   Windows are returned as column vectors.
            %
            %   Usage:
            %      win = DGTtool.getWindow(windowLength,windowName)
            %      [win,diffWin] = DGTtool.getWindow(___)
            %
            %   Acceptable window name (shortcut):
            %      'Hann'           ('h')
            %      'Blackman'       ('b')
            %      '3termC1Nuttall' ('3termC1')
            %      '3termC3Nuttall' ('3termC3')
            %      '4termC1Nuttall' ('4termC1')
            %      '4termC3Nuttall' ('4termC3')
            %      '4termC5Nuttall' ('4termC5')
            %      'Gauss'          ('g')
            %      'Slepian'        ('s')
            %      'Chebyshev'      ('c')
            %
            %   For 'Gauss', 'Slepian' and 'Chebyshev', options can be set.
            %      win = DGTtool.getWindow(___,Name,Value)
            %
            %   Options (width parameters):
            %      'Gauss'     (number, default = 0.02)
            %      'Slepian'   (number, default = 12)
            %      'Chebyshev' (number, default = 340 [dB])
            %
            %   see also DGTtool, DGTtool/DGTtool, setWindow
% [解説 L937] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L938] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                windowLength (1,1) {mustBePositive,mustBeInteger}
% [解説 L939] 処理途中で使う変数 windowName を 任意shape、型 char の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
                windowName   char
% [解説 L940] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                options.Gauss     (1,1) {mustBePositive} = 0.02
% [解説 L941] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                options.Slepian   (1,1) {mustBePositive} = 12
% [解説 L942] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                options.Chebyshev (1,1) {mustBePositive} = 340
            end
            
% [解説 L945] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 winList」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
            winList = DGTtool.windowList;
% [解説 L946] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 winName」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
            winName = validatestring(windowName,winList);
            
% [解説 L948] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 isodd」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            isodd = mod(windowLength,2); % 1 (if odd) or 0 (if even)
% [解説 L949] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 K」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            K = windowLength + isodd;    % always even
            
% [解説 L951] 条件「ismember(winName,winList(1:8))」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ismember(winName,winList(1:8)) % cosine window case
                
% [解説 L953] 値 winName に応じて、入力形式・音源モデル・評価モードなどの処理経路を切り替えます。
                switch winName
% [解説 L954] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                    case 'ForKitamuraHamming'
% [解説 L955] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 c」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                        c = [25; 21]./46;
% [解説 L956] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                    case 'Hann'
% [解説 L957] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 c」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                        c = [0.5; 0.5];
% [解説 L958] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                    case 'Blackman'
% [解説 L959] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 c」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                        c = [0.42; 0.5; 0.08];
% [解説 L960] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                    case '3termC1Nuttall'
% [解説 L961] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 c」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                        c = [0.40897; 0.5; 0.09103];
% [解説 L962] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                    case '3termC3Nuttall'
% [解説 L963] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 c」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                        c = [3; 4; 1]/8;
% [解説 L964] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                    case '4termC1Nuttall'
% [解説 L965] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 c」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                        c = [0.355768; 0.487396; 0.144232; 0.012604];
% [解説 L966] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                    case '4termC3Nuttall'
% [解説 L967] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 c」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                        c = [0.338946; 0.481973; 0.161054; 0.018027];
% [解説 L968] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                    case '4termC5Nuttall'
% [解説 L969] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 c」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                        c = [10; 15; 6; 1]/32;
                end
                
% [解説 L972] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                t = (-(K/2-1):(K/2-1))' / K; % always odd (-0.5 < t < 0.5)
% [解説 L973] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 win」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
                win = cos(2*pi*t.*(0:length(c)-1)) * c;
                
% [解説 L975] 条件「nargout == 2」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if nargout == 2
% [解説 L976] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 diffC」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                    diffC = (0:length(c)-1)' .* c;
% [解説 L977] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 diffWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                    diffWin = -sin(2*pi*t.*(0:length(c)-1)) * diffC;
                end
                
% [解説 L980] 直前までの条件に該当しなかった場合の代替処理へ進みます。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
            else % other case
                
% [解説 L953] 値 winName に応じて、入力形式・音源モデル・評価モードなどの処理経路を切り替えます。
                switch winName
% [解説 L983] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                    case 'Gauss'
% [解説 L984] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                        t = (-(K/2-1):(K/2-1))' / K;
% [解説 L985] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 win」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
                        win = exp(-pi*t.^2/options.Gauss);
% [解説 L986] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                    case 'Slepian'
% [解説 L987] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 win」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
                        win = dpss(K-1,options.Slepian,1);    % Signal Processing Toolbox
% [解説 L988] 候補の最大値を取り、ピークや最良候補を求めます。結果を「処理途中で使う変数 win」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
                        win = win / max(win);
% [解説 L989] 指定された選択値に対応する方式を実行し、同じ関数内で複数のアルゴリズムや設定を切り替えます。
                    case 'Chebyshev'
% [解説 L990] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 win」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
                        win = chebwin(K-1,options.Chebyshev); % Signal Processing Toolbox
                end
                
% [解説 L993] 条件「nargout == 2」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if nargout == 2
% [解説 L994] 条件「isequal(winName,'Gauss')」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                    if isequal(winName,'Gauss')
% [解説 L995] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 diffWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                        diffWin = -exp(-pi*t.^2/options.Gauss).*t/options.Gauss;
% [解説 L996] この処理を実行して現在の計算状態を更新し、次の処理が必要とする状態へ進めます。
                    elseif win(1) < eps * length(win)
% [解説 L997] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 diffWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                        diffWin = DGTtool.computeNumericalDiffWin(win,2^nextpow2(2*length(win)));
% [解説 L998] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 diffWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                        diffWin = diffWin(1:length(win));
% [解説 L999] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
                    else
% [解説 L1000] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 diffWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                        diffWin = [];
                    end
                end
                
            end
            
% [解説 L1006] 条件「~isodd」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isodd % even
% [解説 L1007] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 win」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
                win = [0; win];
% [解説 L1008] 条件「exist('diffWin','var') && ~isempty(diffWin)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if exist('diffWin','var') && ~isempty(diffWin)
% [解説 L1009] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 diffWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                    diffWin = [0; diffWin];
                end
            end
        end
        
% [解説 L1014] この関数は入力を受け取り、dualWin を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function dualWin = computeCanonicalDual(win,shift,FFTnum,sigLen)
            %computeCanonicalDual: Compute canonical dual of given window.
            %   Canonical dual ensures perfect reconstruction of signal.
            %   Signal length must be specified whenever FFTnum < winLen.
            %
            %   Usage:
            %      dualWin = DGTtool.computeCanonicalDual(win,shift,FFTnum)
            %      dualWin = DGTtool.computeCanonicalDual(win,shift,FFTnum,sigLen)
            %
            %   See also pinv
% [解説 L1024] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L1025] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                win    (:,1) {mustBeNumeric}
% [解説 L1026] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                shift  (1,1) {mustBePositive,mustBeInteger}
% [解説 L1027] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                FFTnum (1,1) {mustBePositive,mustBeInteger}
% [解説 L1028] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                sigLen (1,1) {mustBeInteger} = -1
            end
            
% [解説 L1031] 条件「FFTnum >= length(win)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if FFTnum >= length(win)
% [解説 L1032] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 dualWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                dualWin = buffer(win,shift);
% [解説 L1033] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 dualWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                dualWin = dualWin ./ sum(abs(dualWin).^2,2);
% [解説 L1034] 要素数を保ったままshapeを組み替えます。結果を「処理途中で使う変数 dualWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                dualWin = reshape(dualWin(1:length(win)),[],1);
% [解説 L1035] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L1036] 条件「sigLen < length(win)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if sigLen < length(win)
% [解説 L1037] この条件では正しい計算ができないため、理由を示して停止します。
% [意図] 前提を満たさない入力で誤った数値結果を作らないよう、処理を停止するためです。
                    error 'sigLen must be specified whenever FFTnum < winLen.'
                end
% [解説 L1039] 時間信号と時間周波数表現の間を変換します。結果を「時間信号のサンプル数」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
                sigLen = DGTtool.sigLenForFactorDGT(sigLen,shift,FFTnum);
% [解説 L1040] この関数の複数の計算結果を受け取り、c、p、q、d、k、l、r、s にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                [c,~,p,q,d,k,l,r,s] = getConstantsForFacAlg(shift,FFTnum,sigLen/shift);
% [解説 L1041] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 idx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
                idx = getWinIdxForFacAlg(k,l,r,s,c,p,q,d);
                
% [解説 L1043] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 win」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
                win = DGTtool.extendSignalByZeroPad(win,sigLen);
% [解説 L1044] 時間方向の信号を周波数成分へ変換します。結果を「処理途中で使う変数 phi」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                phi = fft(win(idx),[],4);
                
% [解説 L1046] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「正解音源または正解音像の時間信号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                S = pagemtimes(phi,'none',phi,'ctranspose');
% [解説 L1047] 処理途中で使う変数 m を 1:size(phi,4) の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各マイクチャネルを個別に処理した後で、マルチチャネル情報を正しく統合するためです。
                for m = 1:size(phi,4)
% [解説 L1048] 処理途中で使う変数 n を 1:size(phi,3) の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
                    for n = 1:size(phi,3)
% [解説 L1049] 逆行列を明示せず線形方程式を解き、数値的に安定な係数を求めます。結果を「処理途中で使う変数 phi(:,:,n,m)」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                        phi(:,:,n,m) = S(:,:,n,m) \ phi(:,:,n,m);
                    end
                end
                
% [解説 L1053] 周波数成分を時間方向の信号へ戻します。結果を「処理途中で使う変数 phi」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                phi = ifft(phi,[],4,'symmetric');
% [解説 L1054] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 dualWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                dualWin = zeros(size(win));
% [解説 L1055] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 dualWin(idx)」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
                dualWin(idx) = phi;
            end
        end
        
% [解説 L1059] この関数は入力を受け取り、tightWin を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function tightWin = computeCanonicalTight(win,shift,FFTnum,sigLen)
            %computeCanonicalTight: Compute canonical tight window of given window.
            %   Using canonical tight window for both analysis and synthesis
            %   results in perfect reconstruction of signal.
            %   Signal length must be specified whenever FFTnum < winLen.
            %
            %   Usage:
            %      tightWin = DGTtool.computeCanonicalTight(win,shift,FFTnum)
            %      tightWin = DGTtool.computeCanonicalTight(win,shift,FFTnum,sigLen)
            %
            %   See also makeWindowTight
% [解説 L1070] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L1071] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                win    (:,1) {mustBeNumeric}
% [解説 L1072] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                shift  (1,1) {mustBePositive,mustBeInteger}
% [解説 L1073] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                FFTnum (1,1) {mustBePositive,mustBeInteger}
% [解説 L1074] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                sigLen (1,1) {mustBeInteger} = -1
            end
            
% [解説 L1077] 条件「FFTnum >= length(win)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if FFTnum >= length(win)
% [解説 L1078] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 tightWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                tightWin = buffer(win,shift);
% [解説 L1079] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 tightWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                tightWin = tightWin ./ sqrt(sum(abs(tightWin).^2,2));
% [解説 L1080] 要素数を保ったままshapeを組み替えます。結果を「処理途中で使う変数 tightWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                tightWin = reshape(tightWin(1:length(win)),[],1);
% [解説 L1081] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L1082] 条件「sigLen < length(win)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if sigLen < length(win)
% [解説 L1083] この条件では正しい計算ができないため、理由を示して停止します。
% [意図] 前提を満たさない入力で誤った数値結果を作らないよう、処理を停止するためです。
                    error 'sigLen (>= winLen) must be specified whenever FFTnum < winLen.'
                end
% [解説 L1085] 時間信号と時間周波数表現の間を変換します。結果を「時間信号のサンプル数」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
                sigLen = DGTtool.sigLenForFactorDGT(sigLen,shift,FFTnum);
% [解説 L1086] この関数の複数の計算結果を受け取り、c、p、q、d、k、l、r、s にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                [c,~,p,q,d,k,l,r,s] = getConstantsForFacAlg(shift,FFTnum,sigLen/shift);
% [解説 L1087] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 idx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
                idx = getWinIdxForFacAlg(k,l,r,s,c,p,q,d);
                
% [解説 L1089] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 win」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
                win = DGTtool.extendSignalByZeroPad(win,sigLen);
% [解説 L1090] 時間方向の信号を周波数成分へ変換します。結果を「処理途中で使う変数 phi」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                phi = fft(win(idx),[],4);
                
% [解説 L1092] 条件「p == 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if p == 1
% [解説 L1093] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 phi」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                    phi = phi ./ sqrt(sum(abs(phi).^2,[1 2]));
% [解説 L1094] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
                else
% [解説 L1095] 処理途中で使う変数 m を 1:size(phi,4) の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各マイクチャネルを個別に処理した後で、マルチチャネル情報を正しく統合するためです。
                    for m = 1:size(phi,4)
% [解説 L1096] 処理途中で使う変数 n を 1:size(phi,3) の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
                        for n = 1:size(phi,3)
% [解説 L1097] この関数の複数の計算結果を受け取り、U、V にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
                            [U,~,V] = svd(phi(:,:,n,m),'econ');
% [解説 L1098] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 phi(:,:,n,m)」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                            phi(:,:,n,m) = U*V';
                        end
                    end
                end
                
% [解説 L1103] 周波数成分を時間方向の信号へ戻します。結果を「処理途中で使う変数 phi」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                phi = ifft(phi,[],4,'symmetric');
% [解説 L1104] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 tightWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                tightWin = zeros(size(win));
% [解説 L1105] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 tightWin(idx)」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
                tightWin(idx) = phi;
            end
        end
        
% [解説 L1109] この関数は入力を受け取り、diffWin を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function diffWin = computeNumericalDiffWin(win,sigLen)
            %computeNumericalDiffWin: Compute numerical differential.
            %   Differential of window is necessary for reassignment.
            %   If sigLen is specified, output is extended to sigLen.
            %
            %   Usage:
            %      diffWin = DGTtool.computeNumericalDiffWin(win)
            %      diffWin = DGTtool.computeNumericalDiffWin(win,sigLen)
            %
            %   See also reassign, plotReassign
% [解説 L1119] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L1120] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                win    (:,1) {mustBeNumeric}
% [解説 L1121] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                sigLen (1,1) {mustBePositive,mustBeInteger} = length(win)
            end
            
% [解説 L1124] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 origWinLen」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
            origWinLen = length(win);
% [解説 L1125] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 win」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
            win = DGTtool.extendSignalByZeroPad(win,sigLen);
            
% [解説 L1127] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 M」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
            M = floor((length(win)-1)/2);
% [解説 L1128] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 fftIdx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
            fftIdx = ifftshift([zeros(mod(length(win)-1,2)),-M:M]); % mod for odd/even cases
% [解説 L1129] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 fftIdx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
            fftIdx = fftIdx(:)*origWinLen/length(win);
            
% [解説 L1131] 時間方向の信号を周波数成分へ変換します。結果を「処理途中で使う変数 diffWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
            diffWin = ifft(1i*fftIdx.*fft(win),'symmetric'); % spectral method
        end
        
% [解説 L1134] この関数は入力を受け取り、tf,reconstErrorBound を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function [tf,reconstErrorBound] = isdual(win,dualWin,shift,FFTnum)
            %ISDUAL: Check whether two windows are dual of each other.
            %   Dual window pair can perfectly reconstruct signal.
            %   This function can return upper bound of (relative) reconstruction error.
            %
            %   Usage:
            %      trueFalse = DGTtool.isdual(win1,win2,shift,FFTnum)
            %      [trueFalse,relReconstErrorBound] = DGTtool.isdual(___)
            %
            %   See also computeCanonicalDual, computeCanonicalTight
% [解説 L1144] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L1145] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                win     (:,1) {mustBeNumeric}
% [解説 L1146] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                dualWin (:,1) {mustBeNumeric}
% [解説 L1147] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                shift   (1,1) {mustBePositive,mustBeInteger}
% [解説 L1148] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                FFTnum  (1,1) {mustBePositive,mustBeInteger}
            end
            
% [解説 L1151] 候補の最大値を取り、ピークや最良候補を求めます。結果を「処理途中で使う変数 maxLen」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
            maxLen  = max(length(win),length(dualWin));
% [解説 L1152] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 win」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
            win     = DGTtool.zeroPadForFactorDGT(win,    shift,FFTnum,maxLen);
% [解説 L1153] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 dualWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
            dualWin = DGTtool.zeroPadForFactorDGT(dualWin,shift,FFTnum,maxLen);
% [解説 L1154] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 idx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
            idx = getIndexForFactorAlg(FFTnum,shift,length(win)/FFTnum,length(win),1);
            
% [解説 L1156] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 WexlerRaz」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
            WexlerRaz = DGT_factorAlg(dualWin,win,FFTnum,shift,idx);
% [解説 L1157] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 WexlerRaz(1,1)」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
            WexlerRaz(1,1) = WexlerRaz(1,1) - shift;
% [解説 L1158] 複素値または符号付き値を振幅・大きさへ変換します。結果を「処理途中で使う変数 absWR」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            absWR = abs(WexlerRaz)/FFTnum;
% [解説 L1159] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 reconstErrorBound」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
            reconstErrorBound = sum(absWR,'all') + sum(absWR(2:end-1+mod(shift,2),:),'all');
            
% [解説 L1161] 候補の最大値を取り、ピークや最良候補を求めます。結果を「処理途中で使う変数 tf」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            tf = reconstErrorBound < eps(max(abs(win))) * length(win);
        end
        
% [解説 L1164] この関数は入力を受け取り、y を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function y = zeroPad(x,zeroNum)
            %zeroPad: Zero-padding by specifying number of added zero.
            %   This function allows multi-channel signal.
            %
            %   Usage:
            %      x = DGTtool.zeroPad(x,zeroNum)
            %
            %   See also extendSignalByZeroPad, zeroPadForFactorDGT
% [解説 L1172] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L1173] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                x       (:,:) {mustBeNumeric,mustBeSkinny}
% [解説 L1174] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                zeroNum (1,1) {mustBeNonnegative,mustBeInteger}
            end
            
% [解説 L1177] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 newLength」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
            newLength = size(x,1) + zeroNum;
% [解説 L1178] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「分離信号の複素時間周波数表現」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
            y = zeros(newLength,size(x,2));
% [解説 L1179] 処理途中で使う変数 n を 1:size(x,2) の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
            for n = 1:size(x,2)
% [解説 L1180] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「分離信号の複素時間周波数表現」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
                y(:,n) = buffer(x(:,n),newLength);
            end
        end
        
% [解説 L1184] この関数は入力を受け取り、y を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function y = extendSignalByZeroPad(x,outputLength)
            %extendSignalByZeroPad: Zero-padding by specifying length of output.
            %   This function allows multi-channel signal.
            %
            %   Usage:
            %      x = DGTtool.extendSignalByZeroPad(x,outputLength)
            %
            %   See also zeroPad, zeroPadForFactorDGT
            
% [解説 L1193] 条件「outputLength < size(x,1)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if outputLength < size(x,1)
% [解説 L1194] この条件では正しい計算ができないため、理由を示して停止します。
% [意図] 前提を満たさない入力で誤った数値結果を作らないよう、処理を停止するためです。
                error 'Must satisfy outputLength >= size(x,1)'
            end
% [解説 L1196] 時間信号と時間周波数表現の間を変換します。結果を「分離信号の複素時間周波数表現」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
            y = DGTtool.zeroPad(x,outputLength-size(x,1));
        end
        
% [解説 L1199] この関数は入力を受け取り、y を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function y = zeroPadForFactorDGT(x,shift,FFTnum,minLen)
            %zeroPadForFactorDGT: Zero-padding for factorization algorithm.
            %   Number of zeros is determined to satisfy L = aN = bM.
            %   Signal can be further extended by specifying lower bound of length.
            %
            %   Usage:
            %      x = DGTtool.zeroPadForFactorDGT(x,shift,FFTnum)
            %      x = DGTtool.zeroPadForFactorDGT(x,shift,FFTnum,minLen)
            %
            %   See also sigLenForFactorDGT, zeroPad, extendSignalByZeroPad
            
            % This zero-padding is necessary for using DGT_factorAlg.
            % The result may be unnecessarily long for DGT_usualAlg.
% [解説 L1212] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L1213] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                x      (:,:) {mustBeNumeric,mustBeSkinny}
% [解説 L1214] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                shift  (1,1) {mustBePositive,mustBeInteger}
% [解説 L1215] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                FFTnum (1,1) {mustBePositive,mustBeInteger}
% [解説 L1216] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                minLen (1,1) {mustBePositive,mustBeInteger} = size(x,1)
            end
            
% [解説 L1219] 条件「size(x,1) > minLen」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if size(x,1) > minLen
% [解説 L1220] この条件では正しい計算ができないため、理由を示して停止します。
% [意図] 前提を満たさない入力で誤った数値結果を作らないよう、処理を停止するためです。
                error 'Must satisfy minLen >= size(x,1)'
            end
            
% [解説 L1223] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 newLength」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
            newLength = DGTtool.sigLenForFactorDGT(minLen,shift,FFTnum);
% [解説 L1224] 時間信号と時間周波数表現の間を変換します。結果を「分離信号の複素時間周波数表現」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
            y = DGTtool.extendSignalByZeroPad(x,newLength);
        end
        
% [解説 L1227] この関数は入力を受け取り、newLength を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function newLength = sigLenForFactorDGT(xLen,shift,FFTnum)
            %sigLenForFactorDGT: Compute signal length necessary for factorization algorithm.
            %   Obtained length L satisfies L = aN = bM.
            %
            %   Usage:
            %      newLen = DGTtool.sigLenForFactorDGT(sigLen,shift,FFTnum)
            %
            %   See also zeroPadForFactorDGT, extendSignalByZeroPad
% [解説 L1235] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L1236] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                xLen   (1,1) {mustBePositive,mustBeInteger}
% [解説 L1237] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                shift  (1,1) {mustBePositive,mustBeInteger}
% [解説 L1238] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                FFTnum (1,1) {mustBePositive,mustBeInteger}
            end
            
% [解説 L1241] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 c」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            c = lcm(shift,FFTnum);
% [解説 L1242] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 newLength」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
            newLength = ceil(xLen/c) * c;
        end
    end
    
% [解説 L1246] DGT、逆変換、描画、窓生成、設定更新を行うメソッド群を定義します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
    methods % set and get (property access methods)
        
% [解説 L1248] プロパティ変更時に関連状態も更新するsetterメソッドを定義します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        function set.shift(obj,shift)
% [解説 L1249] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L1250] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                obj
% [解説 L1251] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                shift (1,1) {mustBePositive,mustBeInteger}
            end
% [解説 L1253] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.shift」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            obj.shift = shift;
% [解説 L1254] checkRedundancy を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
            checkRedundancy(obj)
% [解説 L1255] checkWindowAndShift を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
            checkWindowAndShift(obj)
% [解説 L1256] deleteInternalParameters を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
            deleteInternalParameters(obj)
        end
        
% [解説 L1259] プロパティ変更時に関連状態も更新するsetterメソッドを定義します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        function set.FFTnum(obj,FFTnum)
% [解説 L1260] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L1261] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                obj
% [解説 L1262] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                FFTnum (1,1) {mustBePositive,mustBeInteger}
            end
% [解説 L1264] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.FFTnum」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
            obj.FFTnum = FFTnum;
% [解説 L1265] checkRedundancy を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
            checkRedundancy(obj)
% [解説 L1266] deleteInternalParameters を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
            deleteInternalParameters(obj)
        end
        
% [解説 L1269] プロパティ変更時に関連状態も更新するsetterメソッドを定義します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        function set.redundancy(obj,redundancy)
% [解説 L1270] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.FFTnum」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
            obj.FFTnum = ceil(obj.shift * redundancy);
% [解説 L1271] deleteInternalParameters を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
            deleteInternalParameters(obj)
        end
        
% [解説 L1274] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        function redundancy = get.redundancy(obj)
% [解説 L1275] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 redundancy」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
            redundancy = obj.FFTnum / obj.shift;
        end
        
% [解説 L1278] プロパティ変更時に関連状態も更新するsetterメソッドを定義します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        function set.win(obj,win)
% [解説 L1279] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L1280] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                obj
% [解説 L1281] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                win {mustBeSkinny}
            end
% [解説 L1283] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.win」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
            obj.win = win;
% [解説 L1284] deleteDualWin を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
            deleteDualWin(obj)
% [解説 L1285] deleteDiffWin を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
            deleteDiffWin(obj)
        end
        
% [解説 L1288] プロパティ変更時に関連状態も更新するsetterメソッドを定義します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
        function set.dualWin(obj,win)
% [解説 L1289] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
            arguments
% [解説 L1290] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                obj
% [解説 L1291] この入力のshapeと許容条件を宣言し、想定外データを計算前に拒否します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                win {mustBeSkinny}
            end
% [解説 L1293] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.dualWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
            obj.dualWin = win;
% [解説 L1294] 条件「~isempty(win)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isempty(win)
% [解説 L1295] checkDualWin を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                checkDualWin(obj)
            end
        end
    end
    
% [解説 L1300] DGT、逆変換、描画、窓生成、設定更新を行うメソッド群を定義します。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
    methods (Hidden) % internal functions
        
% [解説 L1302] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function computeIndexForFactorAlg(obj,segNum,sigLen,chNum)
% [解説 L1303] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.factorIdx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
            obj.factorIdx = getIndexForFactorAlg(obj.shift,obj.FFTnum,segNum,sigLen,chNum);
% [解説 L1304] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
            disp 'index computed (factor DGT)'
        end
        
% [解説 L1307] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function computeIndexForOLA(obj,X,win)
% [解説 L1308] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.OLAindex」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
            obj.OLAindex = getIndexForOLA(X,win,obj.shift);
% [解説 L1309] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
            disp 'index computed (OLA)'
        end
        
% [解説 L1312] この関数は入力を受け取り、x を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function x = zeroPadForFactorAlg(obj,x,newLength)
% [解説 L1313] 時間信号と時間周波数表現の間を変換します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
            x = DGTtool.zeroPadForFactorDGT(x,obj.shift,obj.FFTnum,newLength);
        end
        
% [解説 L1316] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function checkRedundancy(obj)
% [解説 L1317] 条件「~isempty(obj.shift) && ~isempty(obj.FFTnum)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isempty(obj.shift) && ~isempty(obj.FFTnum)
% [解説 L1318] 条件「obj.redundancy <= 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if obj.redundancy <= 1
% [解説 L1319] 不正な入力のまま誤った結果を作らないよう、理由を示して処理を停止します。
                    error(['FFTnum > shift is required for reconstruction! ' ...
                        '(FFTnum = ' num2str(obj.FFTnum) ', shift = ' num2str(obj.shift) ')'])
                end
            end
        end
        
% [解説 L1325] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function checkWindowAndShift(obj)
% [解説 L1326] 条件「~isempty(obj.shift) && ~isempty(obj.win)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isempty(obj.shift) && ~isempty(obj.win)
% [解説 L1327] 条件「obj.shift >= nnz(obj.win)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if obj.shift >= nnz(obj.win)
% [解説 L1328] 不正な入力のまま誤った結果を作らないよう、理由を示して処理を停止します。
                    error(['winLen > shift is required for reconstruction! ' ...
                        '(winLen = ' num2str(nnz(obj.win)) ', shift = ' num2str(obj.shift) ')'])
                end
            end
        end
        
% [解説 L1334] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function deleteInternalParameters(obj)
% [解説 L1335] 条件「~isempty(obj.dualWin)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isempty(obj.dualWin)
% [解説 L1336] deleteDualWin を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                deleteDualWin(obj)
            end
% [解説 L1338] 条件「~isempty(obj.sigLen)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isempty(obj.sigLen)
% [解説 L1339] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.sigLen」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
                obj.sigLen = [];
            end
% [解説 L1341] 条件「~isempty(obj.OLAindex)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isempty(obj.OLAindex)
% [解説 L1342] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.OLAindex」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
                obj.OLAindex = [];
            end
% [解説 L1344] 条件「~isempty(obj.factorIdx)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isempty(obj.factorIdx)
% [解説 L1345] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.factorIdx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
                obj.factorIdx = [];
            end
% [解説 L1347] 条件「~isempty(obj.defConverter)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isempty(obj.defConverter)
% [解説 L1348] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.defConverter」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                obj.defConverter = [];
            end
% [解説 L1350] 条件「~isempty(obj.zeroPhaseConverter)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isempty(obj.zeroPhaseConverter)
% [解説 L1351] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.zeroPhaseConverter」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                obj.zeroPhaseConverter = [];
            end
        end
        
% [解説 L1355] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function deleteDualWin(obj)
% [解説 L1356] 条件「~isempty(obj.dualWin) && obj.isNotWinCalcInDGT」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isempty(obj.dualWin) && obj.isNotWinCalcInDGT
% [解説 L1357] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.dualWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                obj.dualWin = [];
% [解説 L1358] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.isDual」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                obj.isDual = false;
% [解説 L1359] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.isCanonical」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                obj.isCanonical = false;
% [解説 L1360] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                disp 'dualWin deleted'
% [解説 L1361] この処理を実行して現在の計算状態を更新し、次の処理が必要とする状態へ進めます。
            elseif ~obj.isNotWinCalcInDGT
% [解説 L1362] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.isNotWinCalcInDGT」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                obj.isNotWinCalcInDGT = true;
            end
        end
        
% [解説 L1366] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function deleteDiffWin(obj)
% [解説 L1367] 条件「~isempty(obj.diffWin)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if ~isempty(obj.diffWin)
% [解説 L1368] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.diffWin」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                obj.diffWin = [];
% [解説 L1369] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                disp 'diffWin deleted'
            end
        end
        
% [解説 L1373] この関数は入力を受け取り、canonicalDual を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
        function canonicalDual = getCanonicalDualWin(obj,sigLen)
% [解説 L1374] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 canonicalDual」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
            canonicalDual = DGTtool.computeCanonicalDual(obj.win,obj.shift,obj.FFTnum,sigLen);
% [解説 L1375] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.isNotCompDual」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
            obj.isNotCompDual = false;
% [解説 L1376] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.isDual」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
            obj.isDual = true;
% [解説 L1377] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.isCanonical」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
            obj.isCanonical = true;
        end
        
% [解説 L1380] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function checkDualWin(obj)
% [解説 L1381] 条件「obj.isNotCompDual」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if obj.isNotCompDual
% [解説 L1382] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 obj.isDual」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                obj.isDual = DGTtool.isdual(obj.win,obj.dualWin,obj.shift,obj.FFTnum);
% [解説 L1383] setCanonicalFlag を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
                setCanonicalFlag(obj)
% [解説 L1384] この行で、時間周波数変換のクラス構造・入力・表示に必要な処理を行います。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
                disp 'dual window checked'
% [解説 L1385] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L1386] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.isNotCompDual」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                obj.isNotCompDual = true;
            end
        end
        
% [解説 L1390] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function setCanonicalFlag(obj)
% [解説 L1391] 条件「obj.isDual」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
            if obj.isDual
% [解説 L1392] 時間信号と時間周波数表現の間を変換します。結果を「処理途中で使う変数 canonicalDual」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
                canonicalDual = DGTtool.computeCanonicalDual(obj.win,obj.shift,obj.FFTnum,length(obj.dualWin));
% [解説 L1393] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 winError」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
                winError = norm(obj.dualWin - canonicalDual);
% [解説 L1394] 条件「winError < eps(max(abs(canonicalDual))) * length(canonicalDual)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
                if winError < eps(max(abs(canonicalDual))) * length(canonicalDual)
% [解説 L1395] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.isCanonical」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                    obj.isCanonical = true;
% [解説 L1396] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
                else
% [解説 L1397] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.isCanonical」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                    obj.isCanonical = false;
                end
% [解説 L1399] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
            else
% [解説 L1400] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.isCanonical」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
                obj.isCanonical = false;
            end
        end
        
% [解説 L1404] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
        function setFlag_winCalcInDGT(obj)
% [解説 L1405] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 obj.isNotWinCalcInDGT」へ代入します。
% [意図] 解析窓と合成窓の条件を整え、DGT後の信号を正しい振幅と時間配置で再構成するためです。
            obj.isNotWinCalcInDGT = false;
        end
    end
end



% ------------------------------------------------------------
% DGT/IDGT
% ------------------------------------------------------------

% [解説 L1416] この関数は入力を受け取り、X,f,t を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function [X,f,t] = DGT_usualAlg(x,win,shift,FFTnum,paddedSiglen)
% paddedSiglen = ceil(size(x,1)/shift)*shift;
% This computation is performed outside for speed.
% paddedSiglen must be an integer multiple of shift (otherwise, error occurs at zeros).

% [解説 L1421] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 winLen」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
winLen = length(win);
% [解説 L1422] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 segNum」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
segNum = paddedSiglen / shift; % must be integer
% [解説 L1423] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 overlap」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
overlap = winLen - shift;
% [解説 L1424] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 rotNum」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
rotNum = overlap - (paddedSiglen - size(x,1));

% [解説 L1426] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 wx」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
wx = zeros(winLen,segNum,size(x,2));
% [解説 L1427] 処理途中で使う変数 n を 1:size(x,2) の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
for n = 1:size(x,2)
% [解説 L1428] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 wx(:,:,n)」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
    wx(:,:,n) = buffer(x(:,n),winLen,overlap,buffer(x(end-rotNum+1:end,n),overlap));
end
% [解説 L1430] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 wx」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
wx = win .* wx;

% [解説 L1432] 条件「winLen <= FFTnum」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if winLen <= FFTnum
% [解説 L1433] 時間方向の信号を周波数成分へ変換します。結果を「処理途中で使う変数 Xfull」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
    Xfull = fft(wx,FFTnum);
% [解説 L1434] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
else
% [解説 L1435] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 Xfull」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
    Xfull = zeros(FFTnum,segNum,size(x,2));
% [解説 L1436] 処理途中で使う変数 n を 1:size(wx,3) の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
    for n = 1:size(wx,3)
% [解説 L1437] 処理途中で使う変数 m を 1:size(wx,2) の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各マイクチャネルを個別に処理した後で、マルチチャネル情報を正しく統合するためです。
        for m = 1:size(wx,2)
% [解説 L1438] 時間方向の信号を周波数成分へ変換します。結果を「処理途中で使う変数 Xfull(:,m,n)」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
            Xfull(:,m,n) = fft(wrapData(wx(:,m,n),FFTnum));
        end
    end
end
% [解説 L1442] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
X = Xfull(1:floor(FFTnum/2)+1,:,:);

% [解説 L1444] 条件「nargout > 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if nargout > 1
% [解説 L1445] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「現在処理している周波数ビン番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
    f = (0:size(X,1)-1)'/FFTnum;
% [解説 L1446] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 cg」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
    cg = centerOfGravity(win) - 1;
% [解説 L1447] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
    t = (0:size(X,2)-1)*shift + cg - overlap;
end
end

% [解説 L1451] この関数は入力を受け取り、x を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function x = invDGT_usualAlg(X,win,FFTnum,OLAidx)
% Over-lap add (OLA) is performed by accumarray and OLAidx.

% [解説 L1454] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 wx」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
wx = ifft([X; zeros(FFTnum-size(X,1),size(X,2),size(X,3))],'symmetric');

% [解説 L1456] 条件「size(wx,1) == length(win)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if size(wx,1) == length(win)
% [解説 L1457] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 wx」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
    wx = win .* wx;
% [解説 L1458] この処理を実行して現在の計算状態を更新し、次の処理が必要とする状態へ進めます。
elseif size(wx,1) > length(win)
% [解説 L1459] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 wx」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
    wx = win .* wx(1:length(win),:,:);
% [解説 L1460] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
else
% [解説 L1461] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 xrep」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
    xrep = repmat(wx,ceil(length(win)/size(wx,1)),1,1);
% [解説 L1462] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 wx」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
    wx = win .* xrep(1:length(win),:,:);
end

% [解説 L1465] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
x = accumarray(OLAidx(:),wx(:)); % OLA
% [解説 L1466] 要素数を保ったままshapeを組み替えます。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
x = reshape(x,[],size(wx,3));
end

% [解説 L1469] この関数は入力を受け取り、OLAidx を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function OLAidx = getIndexForOLA(X,win,shift)
% [解説 L1470] 信号またはベクトルの長さを取得します。結果を「処理途中で使う変数 winLen」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
winLen = length(win);
% [解説 L1471] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 overlap」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
overlap = winLen - shift;

% [解説 L1473] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 signalLength」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
signalLength = size(X,2) * shift;
% [解説 L1474] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 idx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
idx = uint32(1:signalLength);
% [解説 L1475] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 OLAindex」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
OLAindex = buffer(idx,winLen,overlap,idx(end-overlap+1:end));
% [解説 L1476] 要素数を保ったままshapeを組み替えます。結果を「処理途中で使う変数 OLAidx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
OLAidx = OLAindex + reshape(uint32(signalLength*(0:size(X,3)-1)),1,1,[]);
end



% ------------------------------------------------------------
% DGT/IDGT (Factorization Algorithm)
% ------------------------------------------------------------

% [解説 L1485] この関数は入力を受け取り、X,f,t を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function [X,f,t] = DGT_factorAlg(x,win,shift,FFTnum,in)
% Peter L. Sondergaard
% Efficient algorithms for the discrete Gabor transform with a long FIR window
% Journal of Fourier Analysis and Applications, vol.18, pp.456-470, 2012.

% [解説 L1490] 条件「isNotInteger(size(x,1)/lcm(shift,FFTnum))」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if isNotInteger(size(x,1)/lcm(shift,FFTnum))
% [解説 L1491] この条件では正しい計算ができないため、理由を示して停止します。
% [意図] 前提を満たさない入力で誤った数値結果を作らないよう、処理を停止するためです。
    error 'Length of input signal L must satisfy L = aN = bM.'
end

% [解説 L1494] 時間方向の信号を周波数成分へ変換します。結果を「処理途中で使う変数 phi」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
phi = fft(win(in.wIdx),[],4);

% [解説 L1496] 時間方向の信号を周波数成分へ変換します。結果を「処理途中で使う変数 psi」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
psi = fft(x(in.xIdx),[],4);
% [解説 L1497] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 C」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
C = pagemtimes(phi,'ctranspose',psi,'none');
% [解説 L1498] 周波数成分を時間方向の信号へ戻します。結果を「処理途中で使う変数 C」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
C = ifft(C,[],4,'symmetric');
% [解説 L1499] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 wx」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
wx = C(in.cIdx);

% [解説 L1501] 条件「size(in.cIdx,2) == 1 && size(in.cIdx,3) == 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if size(in.cIdx,2) == 1 && size(in.cIdx,3) == 1
% [解説 L1502] 時間方向の信号を周波数成分へ変換します。結果を「処理途中で使う変数 Xfull」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
    Xfull = fft(wx(:));
% [解説 L1503] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
    X = Xfull(1:floor(FFTnum/2)+1);
% [解説 L1504] ここからは直前の条件が成立しなかった場合の処理で、別の設定・データ状態に対応します。
else
% [解説 L1505] 時間方向の信号を周波数成分へ変換します。結果を「処理途中で使う変数 Xfull」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
    Xfull = fft(wx);
% [解説 L1506] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
    X = Xfull(1:floor(FFTnum/2)+1,:,:);
end

% [解説 L1509] 条件「nargout > 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if nargout > 1
% [解説 L1510] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「現在処理している周波数ビン番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
    f = (0:size(X,1)-1)'/FFTnum;
% [解説 L1511] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 cg」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
    cg = centerOfGravity(win) - 1;
% [解説 L1512] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
    t = mod((0:size(X,2)-1)*shift + cg, size(x,1));
end
end

% [解説 L1516] この関数は入力を受け取り、x を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function x = invDGT_factorAlg(X,win,shift,FFTnum,in)
% Peter L. Sondergaard
% Efficient algorithms for the discrete Gabor transform with a long FIR window
% Journal of Fourier Analysis and Applications, vol.18, pp.456-470, 2012.

% [解説 L1521] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 chNum」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
chNum = size(X,3);
% [解説 L1522] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「時間信号のサンプル数」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
sigLen = size(X,2) * shift;

% [解説 L1524] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
x = zeros(sigLen,size(X,3));
% [解説 L1525] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 C」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
C = zeros([in.q in.q*chNum in.c in.d]);

% [解説 L1527] 時間方向の信号を周波数成分へ変換します。結果を「処理途中で使う変数 phi」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
phi = fft(win(in.wIdx),[],4);

% [解説 L1529] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 wx」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
wx = ifft([X; zeros(FFTnum-size(X,1),size(X,2),size(X,3))],'symmetric');
% [解説 L1530] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 C(in.cIdx)」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
C(in.cIdx) = wx;
% [解説 L1531] 時間方向の信号を周波数成分へ変換します。結果を「処理途中で使う変数 C」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
C = fft(C,[],4);

% [解説 L1533] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 psi」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
psi = pagemtimes(phi,C);
% [解説 L1534] 周波数成分を時間方向の信号へ戻します。結果を「処理途中で使う変数 psi」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
psi = ifft(psi,[],4,'symmetric');

% [解説 L1536] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「複数音源が混ざった観測時間信号」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
x(in.xIdx) = psi;
end

% [解説 L1539] この関数は入力を受け取り、factorIdx を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function factorIdx = getIndexForFactorAlg(shift,FFTnum,segNum,sigLen,chNum)
% [解説 L1540] この関数の複数の計算結果を受け取り、c、ha、p、q、d、k、l、r、s にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[c,ha,p,q,d,k,l,r,s] = getConstantsForFacAlg(shift,FFTnum,segNum);
% [解説 L1541] この関数の複数の計算結果を受け取り、wIdx、xIdx、cIdx にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[wIdx,xIdx,cIdx] = getIndicesForFacAlg(k,l,r,s,c,p,q,d,ha,sigLen,shift,FFTnum,segNum,chNum);
% [解説 L1542] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 factorIdx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
factorIdx = struct('c',c,'d',d,'p',p,'q',q,'wIdx',wIdx,'xIdx',xIdx,'cIdx',cIdx);
end

% [解説 L1545] この関数は入力を受け取り、c,ha,p,q,d,k,l,r,s を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function [c,ha,p,q,d,k,l,r,s] = getConstantsForFacAlg(shift,FFTnum,segNum)
% [解説 L1546] この関数の複数の計算結果を受け取り、c、ha にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[c,ha,~] = gcd(shift,FFTnum);
% [解説 L1547] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 c」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
c = int32(c);

% [解説 L1549] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 p」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
p = shift/c;
% [解説 L1550] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「観測の共分散を単位行列へ近づける白色化行列」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
q = FFTnum/c;
% [解説 L1551] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 d」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
d = segNum/q;

% [解説 L1553] 要素数を保ったままshapeを組み替えます。結果を「処理途中で使う変数 k」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
k = reshape(0:p-1, [],1,1,1);
% [解説 L1554] 要素数を保ったままshapeを組み替えます。結果を「処理途中で使う変数 l」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
l = reshape(0:q-1, 1,[],1,1);
% [解説 L1555] 要素数を保ったままshapeを組み替えます。結果を「AuxFDICA更新や射影に使う共分散・相関行列」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
r = reshape(1:c  , 1,1,[],1);
% [解説 L1556] 要素数を保ったままshapeを組み替えます。結果を「正解音源または正解音像の時間信号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
s = reshape(0:d-1, 1,1,1,[]);
end

% [解説 L1559] この関数は入力を受け取り、wIdx,xIdx,cIdx を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function [wIdx,xIdx,cIdx] = getIndicesForFacAlg(k,l,r,s,c,p,q,d,ha,L,shift,FFTnum,segNum,chNum)
% [解説 L1560] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 wIdx」へ代入します。
% [意図] 周波数ごとの音源分離を表す中心パラメータを初期化または更新するためです。
wIdx = getWinIdxForFacAlg(k,l,r,s,c,p,q,d);
% [解説 L1561] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 xIdx」へ代入します。
% [意図] 周波数ごとに分離行列を適用できる表現へ観測信号を移すためです。
xIdx = getSigIdxForFacAlg(k,l,r,s,p,q,ha,L,shift,FFTnum,chNum);
% [解説 L1562] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 cIdx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
cIdx = getSpecIdxForFacAlg(c,q,d,ha,FFTnum,segNum,chNum);
end

% [解説 L1565] この関数は入力を受け取り、idx を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function idx = getWinIdxForFacAlg(k,l,r,s,c,p,q,d)
% [解説 L1566] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 idx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
idx = r + c*mod(k*q - l*p + s*(p*q), d*p*q);
end

% [解説 L1569] この関数は入力を受け取り、idx を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function idx = getSigIdxForFacAlg(k,l,r,s,p,q,ha,L,shift,FFTnum,chNum)
% [解説 L1570] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 idx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
idx = r + mod(k*FFTnum + s*(p*FFTnum) + l*(ha*shift), L);
% [解説 L1571] 条件「chNum > 1」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if chNum > 1
% [解説 L1572] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 idx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
    idx = repmat(idx,1,chNum) + int32(repelem(L*(0:chNum-1),q));
end
end

% [解説 L1576] この関数は入力を受け取り、idx を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function idx = getSpecIdxForFacAlg(c,q,d,ha,FFTnum,segNum,chNum)
% [解説 L1577] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 sizeC」へ代入します。
% [意図] 実データの軸長をループ範囲と配列shapeへ反映し、次元不一致を防ぐためです。
sizeC = [q q*chNum c d];
% [解説 L1578] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「処理途中で使う変数 idx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
idx = reshape(permute(reshape(int32(1:prod(sizeC)),sizeC),[3 1 4 2]),c,q*d,q,chNum);
% [解説 L1579] 処理途中で使う変数 n を 1:q の範囲で変えながら、同じ処理を各要素へ適用します。
% [意図] 各音源候補を個別に更新・比較し、音源ごとの結果を正しく保持するためです。
for n = 1:q
% [解説 L1580] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 idx(:,:,n,:)」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
    idx(:,:,n,:) = circshift(idx(:,:,n,:),(n-1)*ha,2);
end
% [解説 L1582] 軸順を並べ替え、行列演算が期待する次元配置へ整えます。結果を「処理途中で使う変数 idx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
idx = reshape(permute(idx,[1 3 2 4]),FFTnum,segNum,chNum);
end



% ------------------------------------------------------------
% Phase manipulation
% ------------------------------------------------------------

% [解説 L1591] この関数は入力を受け取り、defConverter を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function defConverter = calculateDefConverter(X,shift,FFTnum)
% [解説 L1592] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「現在処理している周波数ビン番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
f = (0:size(X,1)-1)';
% [解説 L1593] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「現在処理している時間フレーム番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
t = (0:size(X,2)-1) * shift;
% [解説 L1594] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 idx」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
idx = mod(f*t,FFTnum);
% [解説 L1595] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 defConverter」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
defConverter = exp(-2i*pi*idx/FFTnum);
end

% [解説 L1598] この関数は入力を受け取り、zeroPhaseConverter を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function zeroPhaseConverter = calculateZeroPhaseConverter(X,FFTnum,win)
% [解説 L1599] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 rotNum」へ代入します。
% [意図] 現在の音源推定に応じた統計量を作り、補助関数法の安定な分離行列更新に使うためです。
rotNum = centerOfGravity(win) - 1;
% [解説 L1600] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「現在処理している周波数ビン番号」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
f = (0:size(X,1)-1)';
% [解説 L1601] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 zeroPhaseConverter」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
zeroPhaseConverter = exp(2i*pi*f*rotNum/FFTnum);
end



% ------------------------------------------------------------
% Helper functions
% ------------------------------------------------------------

% [解説 L1610] この関数は入力を受け取り、y を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function y = wrapData(x,M)
%wrapData: Simplified version of datawrap in Signal Processing Toolbox.
%   Validation of input arguments is skipped for speed.
%   x must be a vector (:,1) or (1,:)
%   M must be a positive integer (1,1)
% [解説 L1615] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「分離信号の複素時間周波数表現」へ代入します。
% [意図] 現在の分離行列で得られる各音源候補を計算し、次の重み更新や出力に使うためです。
y = sum(buffer(x,M),2);
end

% [解説 L1618] この関数は入力を受け取り、cg を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function cg = centerOfGravity(x)
% [解説 L1619] ここから入力引数の型・shape・既定値・許容条件を宣言します。
% [意図] 不適切な実験条件を計算開始前に検出し、結果の再現性を守るための確認です。
arguments
% [解説 L1620] 複数音源が混ざった観測時間信号 を (:,1)、型 double の入力として受け取ります。
% [意図] この制約により、後続の行列演算へ想定外の次元や型が入るのを防ぎます。
    x (:,1) double
end
% [解説 L1622] 指定した軸を足し合わせ、混合信号または集約統計量を作ります。結果を「処理途中で使う変数 cg」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
cg = sum((1:length(x))'.*x) / sum(x);
end

% [解説 L1625] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
function windowStylePlot(win)
% [解説 L1626] 信号・スペクトログラム・コスト推移を図にし、分離状態や収束挙動を目で確認できるようにします。
plot(win,'linewidth',2)
% [解説 L1627] xticks を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
xticks([])
% [解説 L1628] yticks を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
yticks([])
% [解説 L1629] yline を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
yline(0)
% [解説 L1630] 図へ補助線または枠線を加え、数値位置を読みやすくします。
% [意図] DGTtoolの公開機能と内部状態を整理し、一貫した条件で変換と逆変換を行うためです。
box on
% [解説 L1631] 条件「~isempty(win)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if ~isempty(win)
% [解説 L1632] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
    xlim([1 length(win)])
% [解説 L1633] 図の軸名・範囲・表示形式を整え、何を比較している図か読み取れるようにします。
    ylim(max(abs(win))*1.4*[-1 1])
end
end

% [解説 L1637] この関数は入力を受け取り、tf を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function tf = isNotInteger(x)
% [解説 L1638] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 tf」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
tf = x ~= round(x);
end

% [解説 L1641] この関数は入力を受け取り、tf を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function tf = sizeMismatch(in,X)
% [解説 L1642] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 cIdxSize」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
cIdxSize = size(in.cIdx);
% [解説 L1643] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 cIdxSize(1)」へ代入します。
% [意図] 処理対象の位置集合を明示し、必要なサンプル・周波数・フレームだけを正しい順序で参照するためです。
cIdxSize(1) = floor(cIdxSize(1)/2) + 1;
% [解説 L1644] 入力配列の各軸長を取得し、データに合う配列やループ範囲を決めます。結果を「処理途中で使う変数 tf」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
tf = ~isequal(cIdxSize,size(X));
end

% [解説 L1647] この関数は入力を受け取り、tf を計算して返します。
% [意図] ファイル全体では「時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。」という役割を担当します。
function tf = factorIdxMismatch(shift,FFTnum,segNum,in)
% [解説 L1648] この関数の複数の計算結果を受け取り、c にそれぞれ保存します。
% [意図] 一つの計算で得られる関連結果を分けて保持し、後続の評価・更新・復元でそれぞれ利用するためです。
[c,~,~] = gcd(shift,FFTnum);
% [解説 L1649] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 c」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
c = int32(c);
% [解説 L1650] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 p」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
p = shift/c;
% [解説 L1651] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「観測の共分散を単位行列へ近づける白色化行列」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
q = FFTnum/c;
% [解説 L1652] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 d」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
d = segNum/q;

% [解説 L1654] 後で値を書き込む配列を必要なshapeで事前確保します。結果を「処理途中で使う変数 checkTF」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
checkTF = zeros(5,1);
% [解説 L1655] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 checkTF(1)」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
checkTF(1) = ~isequal(in.c,c);
% [解説 L1656] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 checkTF(2)」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
checkTF(2) = ~isequal(in.d,d);
% [解説 L1657] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 checkTF(3)」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
checkTF(3) = ~isequal(in.p,p);
% [解説 L1658] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 checkTF(4)」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
checkTF(4) = ~isequal(in.q,q);
% [解説 L1659] 右辺の計算結果を、後続処理で再利用できる形に保持します。結果を「処理途中で使う変数 tf」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
tf = any(checkTF);
end

% [解説 L1662] この関数は引数で受け取ったデータを処理します。
% [意図] 時間信号と時間周波数表現を相互変換するDGT/STFT系ツールです。
function mustBeSkinny(x)
% The first dimension of input matrix x must be greater than the second dimension.
% [解説 L1664] 条件「size(x,1) ~= 0」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
if size(x,1) ~= 0
% [解説 L1665] 条件「size(x,1) < size(x,2)」を確認し、成立した場合だけ直後の処理を行います。
% [意図] 設定またはデータ状態に合った計算経路を選び、無効な処理を避けます。
    if size(x,1) < size(x,2)
% [解説 L1666] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 eidType」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
        eidType = 'mustBeSkinny:notSkinny';
% [解説 L1667] 複素共役を含む積を作り、複素信号の内積・共分散を正しく計算します。結果を「処理途中で使う変数 msgType」へ代入します。
% [意図] DGT/STFTと逆変換の次段階で再利用し、窓・信号長・軸配置の整合を保つためです。
        msgType = 'Input matrix must be skinny (tall), i.e., size(x,1) >= size(x,2)';
% [解説 L1668] throwAsCaller を呼び出し、この処理段階で必要な計算または表示を実行します。
% [意図] 入力shape、窓長、シフト量、FFT数の整合を確認し、DGTと逆変換が対応する条件を保つためです。
        throwAsCaller(MException(eidType,msgType))
    end
end
end
