clear
clc
close all

addpath('.\\online-calculators-v3\\')
addpath('.\\Matlab MCMC ensemble sampler\\')
addpath('.\\CosmoTools\\')

%% USER INPUT ----------------------------------------------------------- %

algorithm = 'hmc';   % Choose algorithm: 'hmc' or 'gwmcmc'
profile = 'quick';   % Choose algorithm profile: 'quick', 'balanced', or 'robust', for manual adjustment check inversion_build_config.m
export = false;
filetag = 'WC_soilmix';    % Use 'test' to run test scenarios
zm = 200;            % soil mixing depth in cm (0 = no mixing)

% files 
DEM  = GRIDobj('.\data\crete_clipped_utm.tif');
file = 'data\WCdata_RFO.xlsx'; % AMS data

% Priors 
T = [1, 6e3];
E_step  = [10, 5e3];
E_spike = [10, 3e2];
LOSS = [1, 250];
CHG = [1, 100];

cfg = inversion_build_config(filetag, algorithm, profile);

% Scenario selection in main script using true/false flags.
allScenarios = {'step', 'samestep', 'samebackground_step', 'samebackground_samestep', ...
    'spike', 'samespike', 'samebackground_spike', 'samebackground_samespike'};
runScenario = [true, false, false, false, true, false, false, false];
cfg.scenarios = allScenarios(runScenario);
cfg.pause = true; % do you want to pause at the plotting stage, before computing the next scenario?

%% Basin geometry and sample metadata
SAMS = cosmosampleread(file);
SAMS = cosmowatersheds(SAMS,DEM);

lat = arrayfun(@(x) median(x.WSLat), SAMS);
lon = arrayfun(@(x) median(x.WSLon), SAMS);
alt = arrayfun(@(x) median(x.WSDEM.Z(:), 'omitnan'), SAMS);

for i = 1:numel(cfg.scenarios)
    scenario = cfg.scenarios{i};
    data = readtable(file);

    %% Observations
    Nobs = [data.N10; data.N14; data.N26];
    dNobs = [data.N10sigma; data.N14sigma; data.N26sigma];

    Nlogical = [~isnan(data.N10) ~isnan(data.N14) ~isnan(data.N26)];

    %% Scenario-specific priors
    if strcmp(scenario, 'step')
        E = E_step;
    else
        E = E_spike;
    end

    [prior_range, var_names] = make_prior_and_varnames( ...
        scenario, T, E, LOSS, CHG, length(data.N10), cfg.nsteps);

    if strcmp(scenario, 'step')
        prior_range(2:11,2) = 300;   % limit erosion in step model erosion rate 1
    end

    %% Constants, production, and forward model
    consts = make_constants();
    sp = sample_parameters(lat, lon, alt, consts);
    forward_model = @(m) Nforward_wrapper(m, sp, consts, zm, scenario, cfg.nsteps, Nlogical);

    %% Initial model
    nWalks = cfg.(cfg.algorithm).nWalks;
    mini = initialmodel_flatprior(prior_range, nWalks);

    %% Likelihood
    lognormpdf = @(x,mu,sigma)-0.5*((x-mu)./sigma).^2 - log(sqrt(2*pi).*sigma);
    logLikeFn = @(m) sum(lognormpdf(Nobs(Nlogical), forward_model(m), dNobs(Nlogical)));

    %% Sample posterior
    tic
    [models, logLikeStore, samplerInfo] = inversion_run_sampler( ...
        cfg.algorithm, prior_range, var_names, mini, logLikeFn, cfg.(cfg.algorithm), scenario, cfg.nsteps);
    runtimeSeconds = toc;

    %% Best-fit model and status
    [best_model, best_model_like] = inversion_select_best_model(models, logLikeStore);
    best_pred = forward_model(double(best_model));

    fprintf('\n%s | %s | %.2fs\n', upper(cfg.algorithm), scenario, runtimeSeconds);
    if strcmp(cfg.algorithm, 'hmc')
        fprintf('Mean HMC acceptance ratio: %.3f\n', mean(samplerInfo.accRatio));
    end
    fprintf('Best model log-likelihood: %.3f\n', best_model_like);

    %% Plots
    h1 = autocorrelationplot(models);
    h2 = chainplot(models, var_names, prior_range);
    h3 = ecornerplot(models, 'ks', true, 'color', [.3 .3 .3], ...
        'name', var_names, 'bestmodel', best_model);
    h4 = barplot_parameters(models, var_names, prior_range, 'bestmodel', best_model);
    h5 = conc_modelledVSobserved(best_pred, data.N10, data.N10sigma, data.N14, data.N14sigma);

    %% Export
    if  export
        outDir = './output/WC_HMC_soilmix';
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end

        base = [outDir '/' cfg.filetag '_' scenario '_' cfg.algorithm];

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
            'best_model', 'best_model_like', 'best_pred', 'prior_range', 'var_names', '-v7.3');

        clear h1 h2 h3 h4 h5
    end

    disp([scenario ' completed'])
    if cfg.pause
        disp('Excecution paused. Press button to continue.')
        pause()
    end
    close all
end
