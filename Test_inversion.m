clear
clc
close all
addpath(genpath(pwd))

%% USER INPUT ----------------------------------------------------------- %

profile = 'quick';   % Choose algorithm profile: 'quick', 'balanced', or 'robust' for manual adjustment check inversion_build_config.m
export = false;
filetag = 'test';    % Use 'test' to run test scenarios

% Priors ------------------------------------------------------------------
T = [1e2, 6e3];
E = [0,   5e3];
LOSS = [0, 200];
CHG = [1, 50];

cfg = inversion_build_config(filetag, profile);

% Scenario selection. Select which scenarios to run:
allScenarios = {'step', 'samestep', 'samebackground_step', 'samebackground_samestep', ...
    'spike', 'samespike', 'samebackground_spike', 'samebackground_samespike',...
    'curve'};
runScenario = [false, false, false, false,...
    false, true, true, true,...
    false];
cfg.scenarios = allScenarios(runScenario);
cfg.pause = true; % do you want to pause at the plotting stage, before computing the next scenario?

% END USER INPUT ---------------------------------------------------------%

% loop through scenarios 
for i = 1:numel(cfg.scenarios)
    scenario = cfg.scenarios{i};

    %% Make test data
    n = cfg.test.nSamples;
    tdata = make_test_data(scenario, n);
    Nlogical = [true(n,2) false(n,1)];

    %% Priors
    [prior_range, var_names] = make_prior_and_varnames( ...
        scenario, T, E, LOSS, CHG, n, tdata.steps);

    %% Constants and production
    consts = make_constants();
    sp = sample_parameters(tdata.lat, tdata.lon, tdata.altitude, consts);

    %% Initial model
    nWalks = cfg.hmc.nWalks;
    mini = initialmodel_flatprior(prior_range, nWalks, tdata.steps);

    %% Forward model and synthetic observations
    forward_model = @(m) Nforward_wrapper(m, sp, consts, scenario, tdata.steps, Nlogical);

    mtest = [tdata.t'; tdata.e'; tdata.changeVariable'];
    testObs = forward_model(mtest);
    sigmaObs = testObs * 0.08;

    %% Likelihood
    lognormpdf = @(x,mu,sigma)-0.5*((x-mu)./sigma).^2 - log(sqrt(2*pi).*sigma);
    logLikeFn = @(m) sum(lognormpdf(testObs, forward_model(m), sigmaObs));

    %% Sample posterior
    tic
    [models, logLikeStore, samplerInfo] = inversion_run_sampler( ...
        prior_range, var_names, mini, logLikeFn, cfg.hmc, scenario, tdata.steps);
    runtimeSeconds = toc;

    %% Best-fit model and quick diagnostics
    [best_model, best_model_like] = inversion_select_best_model(models, logLikeStore);
    best_pred = forward_model(double(best_model));

    true_model_like = logLikeFn(mtest);
    fprintf('\nHMC | %s | %.2fs\n', scenario, runtimeSeconds);
    fprintf('Mean HMC acceptance ratio: %.3f\n', mean(samplerInfo.accRatio));
    fprintf('Best model log-likelihood: %.3f\n', best_model_like);
    fprintf('True model log-likelihood: %.3f\n', true_model_like);
    fprintf('Log-likelihood gap (true-best): %.3f\n', true_model_like - best_model_like);

    %% Plots

    h1 = autocorrelationplot(models);
    h2 = chainplot(models, var_names, prior_range, mtest);
    h3 = ecornerplot(models, 'ks', true, 'color', [.3 .3 .3], ...
        'name', var_names, 'bestmodel', best_model, 'truevals', mtest);
    h4 = barplot_parameters(models, var_names, prior_range, 'bestmodel', best_model, 'truevals', mtest);
    h5 = conc_modelledVSobserved(best_pred, testObs(1:n), sigmaObs(1:n), testObs(n+1:end), sigmaObs(n+1:end));

    %% Export
    if export
        base = ['./output/' cfg.filetag '_' scenario];
        exportgraphics(h1, [base '_autocorrelation.png'], 'Resolution', 300)
        exportgraphics(h2, [base '_chains.png'], 'Resolution', 300)
        exportgraphics(h3, [base '_cornerplot.png'], 'Resolution', 300)
        exportgraphics(h4, [base '_barplot.png'], 'Resolution', 300)
        exportgraphics(h5, [base '_datafit.png'], 'Resolution', 300)

        exportgraphics(h1, [base '_autocorrelation.pdf'], 'ContentType', 'vector')
        exportgraphics(h2, [base '_chains.pdf'], 'ContentType', 'vector')
        exportgraphics(h3, [base '_cornerplot.pdf'], 'ContentType', 'vector')
        exportgraphics(h4, [base '_barplot.pdf'], 'ContentType', 'vector')
        exportgraphics(h5, [base '_datafit.pdf'], 'ContentType', 'vector')

        save([base '_workspace.mat'], ...
            'cfg', 'scenario', 'models', 'logLikeStore', 'samplerInfo', ...
            'best_model', 'best_model_like', 'best_pred', 'mtest', ...
            'testObs', 'sigmaObs', 'prior_range', 'var_names', '-v7.3');

        clear h1 h2 h3 h4 h5
    end

    disp([scenario ' completed'])

    if cfg.pause
        disp('Excecution paused. Press button to continue.')
        pause()
    end
    close all
end
