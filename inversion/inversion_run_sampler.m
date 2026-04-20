function [models, logLikeStore, samplerInfo] = inversion_run_sampler(prior_range, var_names, mini, logLikeFn, samplerCfg, scenario, nsteps)
% Run HMC and return consistent output arrays.

lower_bounds = prior_range(:,1);
upper_bounds = prior_range(:,2);

if samplerCfg.useTargetedStarts
    % This is only important for samestep: It makes sure change values
    % above and below 1 work. I noticed that some chains can get stuck in
    % change values < 1 and not make it over the likelihood ridge. The
    % targeted start makes sure both side of this posterior ridge (change=1) get
    % explored. Could be removed...Rarely important
    mini = apply_targeted_multistart(mini, prior_range, scenario, nsteps);
end

% Positive lower-bounded parameters are sampled in log space by default.
logScaleMask = lower_bounds > 0;

nWalks = samplerCfg.nWalks;
numSamples = samplerCfg.numSamples;
nParams = size(prior_range, 1);

models = nan(nParams, nWalks, numSamples);
logLikeStore = -inf(2, nWalks, numSamples);
accRatio = nan(nWalks, 1);
endPoints = nan(nParams, nWalks);
hmcSamplers = cell(nWalks, 1);
tuningInfo = cell(nWalks, 1);

stepSize = samplerCfg.stepSizeDefault;
if isfield(samplerCfg, 'stepSizeByScenario') && isfield(samplerCfg.stepSizeByScenario, scenario)
    stepSize = samplerCfg.stepSizeByScenario.(scenario);
end

for wix = 1:nWalks
    startModel = mini(:,wix);
    startPoint = bounded_to_unbounded(startModel, lower_bounds, upper_bounds, logScaleMask);
    logpdf = @(u) transformed_logpdf(u, lower_bounds, upper_bounds, logScaleMask, logLikeFn);

hmc = hmcSampler(logpdf, startPoint, ...
    'UseNumericalGradient', true, ...
    'CheckGradient', false, ...
    'VariableNames', var_names, ...
    'StepSize', stepSize, ...
    'NumSteps', samplerCfg.numLeapfrogSteps, ...
    'StepSizeTuningMethod', 'dual-averaging', ...
    'MassVectorTuningMethod', 'hessian');

    if samplerCfg.useTuning
        [hmcSamplers{wix}, tuningInfo{wix}] = tuneSampler(hmc, ...
            'NumStepSizeTuningIterations', samplerCfg.numTune, ...
            'VerbosityLevel', 0);
    else
        hmcSamplers{wix} = hmc;
        tuningInfo{wix} = struct();
    end

    [chainU, endpointU, accRatio(wix)] = drawSamples(hmcSamplers{wix}, ...
        'Burnin', samplerCfg.numBurnin, ...
        'NumSamples', samplerCfg.numSamples, ...
        'ThinSize', samplerCfg.thinSize, ...
        'VerbosityLevel', 0);

    chainM = bounded_from_unbounded(chainU', lower_bounds, upper_bounds, logScaleMask);
    models(:,wix,:) = reshape(chainM, size(chainM,1), 1, size(chainM,2));
    endPoints(:,wix) = bounded_from_unbounded(endpointU, lower_bounds, upper_bounds, logScaleMask);

    for six = 1:samplerCfg.logEvalStride:numSamples
        logLikeStore(1,wix,six) = 0;
        logLikeStore(2,wix,six) = logLikeFn(models(:,wix,six));
    end
end

models = single(models);
logLikeStore = single(logLikeStore);
samplerInfo = struct();
samplerInfo.algorithm = 'hmc';
samplerInfo.accRatio = accRatio;
samplerInfo.endPoints = endPoints;
samplerInfo.hmcSamplers = hmcSamplers;
samplerInfo.tuningInfo = tuningInfo;
samplerInfo.logScaleMask = logScaleMask;
samplerInfo.stepSize = stepSize;

end

%% helper functions

% swith between unbounded lod-space and bounded linear space
function u = bounded_to_unbounded(m, lower_bounds, upper_bounds, logScaleMask)
u = zeros(size(m));
for ix = 1:size(m,1)
    if logScaleMask(ix)
        lo = log(lower_bounds(ix));
        hi = log(upper_bounds(ix));
        z = (log(m(ix,:)) - lo) / (hi - lo);
    else
        z = (m(ix,:) - lower_bounds(ix)) / (upper_bounds(ix) - lower_bounds(ix));
    end
    z = min(max(z, 1e-10), 1 - 1e-10);
    u(ix,:) = log(z ./ (1 - z));
end
end


function m = bounded_from_unbounded(u, lower_bounds, upper_bounds, logScaleMask)
s = 1 ./ (1 + exp(-u));
m = zeros(size(u));
for ix = 1:size(u,1)
    if logScaleMask(ix)
        lo = log(lower_bounds(ix));
        hi = log(upper_bounds(ix));
        m(ix,:) = exp(lo + (hi - lo) * s(ix,:));
    else
        m(ix,:) = lower_bounds(ix) + (upper_bounds(ix) - lower_bounds(ix)) * s(ix,:);
    end
end
end


function lpdf = transformed_logpdf(u, lower_bounds, upper_bounds, logScaleMask, logLikeFn)
m = bounded_from_unbounded(u, lower_bounds, upper_bounds, logScaleMask);

% Stable evaluation of log(sigmoid(u)) and log(1-sigmoid(u)) without a
% separate softplus helper.
log_sigmoid = -(log1p(exp(-abs(-u))) + max(-u,0));
log_one_minus_sigmoid = -(log1p(exp(-abs(u))) + max(u,0));

log_jacobian_terms = zeros(size(u));
for ix = 1:size(u,1)
    if logScaleMask(ix)
        lo = log(lower_bounds(ix));
        hi = log(upper_bounds(ix));
        log_jacobian_terms(ix,:) = log(m(ix,:)) + log(hi - lo) + log_sigmoid(ix,:) + log_one_minus_sigmoid(ix,:);
    else
        span = upper_bounds(ix) - lower_bounds(ix);
        log_jacobian_terms(ix,:) = log(span) + log_sigmoid(ix,:) + log_one_minus_sigmoid(ix,:);
    end
end

lpdf = logLikeFn(m) + sum(log_jacobian_terms, 'all');
end


function mini = apply_targeted_multistart(mini, prior_range, scenario, nsteps)
if ~contains(scenario, 'samestep')
    return
end

chg_idx = (size(prior_range,1)-nsteps+1):size(prior_range,1);
chg_lo = prior_range(chg_idx,1);
chg_hi = prior_range(chg_idx,2);

for wix = 1:size(mini,2)
    frac = (wix - 1) / max(1, size(mini,2)-1);
    if frac < 0.5
        lo = chg_lo;
        hi = chg_lo + 0.2 * (chg_hi - chg_lo);
    else
        lo = chg_lo + 0.6 * (chg_hi - chg_lo);
        hi = chg_hi;
    end
    mini(chg_idx,wix) = lo + rand(size(chg_lo)) .* (hi - lo);
end
end
