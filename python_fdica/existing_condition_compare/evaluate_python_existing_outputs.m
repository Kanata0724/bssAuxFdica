function evaluate_python_existing_outputs()
% Evaluate Python estimates with existing MATLAB bss_eval_sources.

repo = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(repo);
addpath(fullfile(repo, "bss_eval"));

outDir = fullfile(repo, "output", "fdica_existing_conditions_compare");
matlabData = load(fullfile(outDir, "matlab_existing_debug.mat"), "refSig", "obsSig", "params");
pythonData = load(fullfile(outDir, "python_existing_debug.mat"));

refSig = matlabData.refSig;
obsSig = matlabData.obsSig;
params = matlabData.params;

[inSdr, inSir, inSar] = bss_eval_sources(repmat(obsSig(:, params.refMic), [1, params.nSrc]).', refSig.');
names = ["current_none", "current_cor", "current_ips", "matlab_like_none", "matlab_like_cor", "matlab_like_ips"];
fields = ["estSigCurrentNone", "estSigCurrentCor", "estSigCurrentIps", "estSigMatlabLikeNone", "estSigMatlabLikeCor", "estSigMatlabLikeIps"];

metrics.input.sdr = inSdr;
metrics.input.sir = inSir;
metrics.input.sar = inSar;
for i = 1:numel(names)
    est = pythonData.(fields(i));
    [sdr, sir, sar, perm] = bss_eval_sources(est.', refSig.');
    metrics.(names(i)).sdr = sdr;
    metrics.(names(i)).sir = sir;
    metrics.(names(i)).sar = sar;
    metrics.(names(i)).perm = perm;
end

save(fullfile(outDir, "python_existing_bss_eval.mat"), "metrics", "-v7");
end
