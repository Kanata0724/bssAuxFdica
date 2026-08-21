function export_matlab_existing_debug()
% Export MATLAB AuxFDICA values using the existing main.m conditions.
% This is a comparison tool and does not modify production implementation.

repo = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(repo);
addpath(fullfile(repo, "bss_eval"));

params.dataNo = 1;
params.seed = 1;
params.fftSize = 4096;
params.shiftSize = params.fftSize / 2;
params.nSrc = 2;
params.nIter = 50;
params.isWhiten = false;
params.srcModel = 'LAP';
params.refMic = 1;
params.permSolver = 'COR';
params.isDraw = true;
params.isFilt = false;
params.micPos = [0, 0.0566];
params.isPowRatio = true;
params.typeCor = 'Gl+Lo';
params.deltaFreq = 3;
params.ratioFreq = 2;

rng(params.seed);
[dirPath, fileName] = getInputFileNames(params.dataNo);
srcSig = [];
for iSrc = 1:params.nSrc
    [sig, fs] = audioread(fullfile(repo, dirPath, fileName(iSrc)));
    srcSig(:, :, iSrc) = sig;
end
obsSig = sum(srcSig, 3);
peakVal = max(abs(obsSig), [], "all");
if peakVal > 1
    obsSig = 0.99 * obsSig / peakVal;
    srcSig = 0.99 * srcSig / peakVal;
end
refSig = squeeze(srcSig(:, params.refMic, :));
sigLen = size(obsSig, 1);

F = DGTtool("windowName", "b", "windowLength", params.fftSize, "windowShift", params.shiftSize);
obsSpec = F.DGT(obsSig);
if params.isWhiten
    obsSpecInput = local_whitening_debug(obsSpec, params.nSrc);
else
    obsSpecInput = obsSpec(:, :, 1:params.nSrc);
end

[estSpecFdica, demixMat, cost, debug] = local_aux_fdica_debug( ...
    obsSpecInput, params.nIter, params.srcModel);
[estSpecFixed, demixMatFixed] = local_projection_back_debug( ...
    estSpecFdica, obsSpec(:, :, params.refMic), demixMat);

estSigNone = F.pinv(estSpecFixed);
estSigNone = estSigNone(1:sigLen, :);

[estSpecCor, permCor] = permSolverCor( ...
    estSpecFixed, params.isPowRatio, params.typeCor, params.deltaFreq, params.ratioFreq);
estSigCor = F.pinv(estSpecCor);
estSigCor = estSigCor(1:sigLen, :);

srcSpec = F.DGT(squeeze(srcSig(:, params.refMic, :)));
[estSpecIps, permIps] = permSolverIps(estSpecFixed, srcSpec);
estSigIps = F.pinv(estSpecIps);
estSigIps = estSigIps(1:sigLen, :);

[inSdr, inSir, inSar] = bss_eval_sources(repmat(obsSig(:, params.refMic), [1, params.nSrc]).', refSig.');
[noneSdr, noneSir, noneSar, nonePerm] = bss_eval_sources(estSigNone.', refSig.');
[corSdr, corSir, corSar, corPerm] = bss_eval_sources(estSigCor.', refSig.');
[ipsSdr, ipsSir, ipsSar, ipsPerm] = bss_eval_sources(estSigIps.', refSig.');

metrics.input.sdr = inSdr;
metrics.input.sir = inSir;
metrics.input.sar = inSar;
metrics.none.sdr = noneSdr;
metrics.none.sir = noneSir;
metrics.none.sar = noneSar;
metrics.none.perm = nonePerm;
metrics.cor.sdr = corSdr;
metrics.cor.sir = corSir;
metrics.cor.sar = corSar;
metrics.cor.perm = corPerm;
metrics.ips.sdr = ipsSdr;
metrics.ips.sir = ipsSir;
metrics.ips.sar = ipsSar;
metrics.ips.perm = ipsPerm;

outDir = fullfile(repo, "output", "fdica_existing_conditions_compare");
if ~isfolder(outDir)
    mkdir(outDir);
end
save(fullfile(outDir, "matlab_existing_debug.mat"), ...
    "params", "fs", "srcSig", "obsSig", "refSig", "obsSpec", "obsSpecInput", ...
    "estSpecFdica", "demixMat", "cost", "debug", ...
    "estSpecFixed", "demixMatFixed", ...
    "estSigNone", "estSigCor", "estSigIps", "permCor", "permIps", ...
    "metrics", "-v7");
end

function Y = local_whitening_debug(X, N)
[I, J, ~] = size(X, [1, 2, 3]);
Y = zeros(I, J, N);
Xp = permute(X, [3, 2, 1]);
for i = 1:I
    Xi = Xp(:, :, i);
    V = Xi * (Xi') / J;
    [P, D] = eig(V);
    [~, idx] = sort(diag(D), "descend");
    D = D(idx, idx);
    P = P(:, idx);
    dP = P(:, 1:N);
    Yi = sqrt(D) \ (dP') * Xi;
    Y(i, :, :) = Yi.';
end
end

function [Y, W, cost, debug] = local_aux_fdica_debug(X, nIter, srcModel)
[I, J, M] = size(X, [1, 2, 3]);
N = M;
E = repmat(eye(M), [1, 1, I]);
W = E;
Y = X;
Xp = permute(X, [3, 2, 1]);
Xph = pagectranspose(Xp);
Yp = permute(Y, [3, 2, 1]);
cost = zeros(nIter + 1, 1);
cost(1, 1) = local_calc_fdica_cost_debug(Yp, W, srcModel, I, J);

debug.Y0 = Y;
debug.W0 = W;
debug.V_first = [];
debug.wn_first = [];
debug.norm_first = [];
debug.Y_after_iter1 = [];
debug.W_after_iter1 = [];
debug.Y_after_iter10 = [];
debug.W_after_iter10 = [];
debug.Y_after_iter25 = [];
debug.W_after_iter25 = [];

for iIter = 1:nIter
    if strcmp(srcModel, 'LAP')
        Rp = max(abs(Yp), 10000 * eps);
    elseif strcmp(srcModel, 'TVG')
        Rp = max(abs(Yp).^2, 10000 * eps);
    else
        error("unsupported source model");
    end
    invRp = 1 ./ Rp;
    for n = 1:N
        D = repmat(invRp(n, :, :), [M, 1, 1]);
        Vk = pagemtimes(D .* Xp, Xph) / J;
        wn = pagemldivide(pagemtimes(W, Vk), E(:, n, :));
        normVal = sqrt(pagemtimes(pagemtimes(wn, "ctranspose", Vk, "none"), wn));
        wn = wn ./ normVal;
        if iIter == 1 && n == 1
            debug.V_first = Vk(:, :, 1);
            debug.wn_first = wn(:, :, 1);
            debug.norm_first = normVal(:, :, 1);
        end
        wnh = pagectranspose(wn);
        Yp(n, :, :) = pagemtimes(wnh, Xp);
        W(n, :, :) = wnh;
    end
    cost(iIter + 1, 1) = local_calc_fdica_cost_debug(Yp, W, srcModel, I, J);
    if iIter == 1
        debug.Y_after_iter1 = permute(Yp, [3, 2, 1]);
        debug.W_after_iter1 = W;
    end
    if iIter == 10
        debug.Y_after_iter10 = permute(Yp, [3, 2, 1]);
        debug.W_after_iter10 = W;
    end
    if iIter == 25
        debug.Y_after_iter25 = permute(Yp, [3, 2, 1]);
        debug.W_after_iter25 = W;
    end
end
Y = permute(Yp, [3, 2, 1]);
debug.Y_final = Y;
debug.W_final = W;
end

function costVal = local_calc_fdica_cost_debug(Yp, W, srcModel, I, J)
detW = zeros(I, 1);
for i = 1:I
    detW(i, 1) = det(W(:, :, i));
end
if strcmp(srcModel, 'LAP')
    costVal = sum(abs(Yp), "all") - 2 * J * sum(log(abs(detW)));
elseif strcmp(srcModel, 'TVG')
    costVal = sum(log(max(abs(Yp).^2, eps)), "all") - 2 * J * sum(log(abs(detW)));
else
    error("unsupported source model");
end
end

function [fixY, fixW] = local_projection_back_debug(Y, S, W)
Yp = permute(Y, [3, 2, 1]);
Sp = permute(S, [3, 2, 1]);
Wp = permute(W, [4, 1, 2, 3]);
Yph = pagectranspose(Yp);
YpYph = pagemtimes(Yp, Yph);
YphOnYpYph = pagemrdivide(Yph, YpYph);
A = pagemtimes(Sp, YphOnYpYph);
Ap = permute(A, [1, 2, 4, 3]);
Ypp = permute(Yp, [4, 1, 2, 3]);
fixY = Ap .* Ypp;
fixY = permute(fixY, [4, 3, 2, 1]);
fixW = Ap .* Wp;
fixW = permute(fixW, [2, 3, 4, 1]);
end
