function cfg = inversion_build_config(filetag, profile)
% Build layered defaults for unified inversion runners.
% filetag = "test" or "yourProjectname"
% profile = 'quick' , 'robust', 'balanced'

cfg = struct();
cfg.filetag = lower(filetag);
cfg.profile = lower(profile);

if strcmp(cfg.filetag, 'test') 
    cfg.nsteps = 1;
    cfg.test.nSamples = 7;
else % if not 'test' scenario
    cfg.nsteps = 1;
end

cfg.hmc.nWalks = 8;
cfg.hmc.useTuning = false;
cfg.hmc.numTune = 20;
cfg.hmc.numBurnin = 5000;
cfg.hmc.numSamples = 20000;
cfg.hmc.thinSize = 1;
cfg.hmc.numLeapfrogSteps = 10;
cfg.hmc.stepSizeDefault = 0.04;
cfg.hmc.stepSizeByScenario = struct( ...
    'step', 0.08, ...
    'samestep', 0.04, ...
    'samebackground_step', 0.04, ...
    'samebackground_samestep', 0.04, ...
    'spike', 0.08, ...
    'samespike', 0.1, ...
    'samebackground_spike', 0.06, ...
    'samebackground_samespike', 0.06);
cfg.hmc.logEvalStride = 1;
cfg.hmc.useLogScalePositive = true;
cfg.hmc.useTargetedStarts = true;

switch cfg.profile
    case 'quick'
        cfg.hmc.nWalks = min(cfg.hmc.nWalks, 2);
        cfg.hmc.useTuning = false;
        cfg.hmc.numTune = 0;
        cfg.hmc.numBurnin = 500;
        cfg.hmc.numSamples = 2000;
        cfg.hmc.thinSize = 1;
        cfg.hmc.numLeapfrogSteps = 6;
        cfg.hmc.logEvalStride = 10;

    case 'robust'
        % cfg.hmc.useTuning = true;
        cfg.hmc.numTune = 80;
        cfg.hmc.numBurnin = 10000;
        cfg.hmc.numSamples = 50000;
        cfg.hmc.logEvalStride = 1;

    case 'balanced'
        % keep defaults

    otherwise
        error('Unknown profile "%s". Use "quick", "balanced", or "robust".', cfg.profile);
end

cfg.filetag = [cfg.filetag '_hmc_' cfg.profile];
end
