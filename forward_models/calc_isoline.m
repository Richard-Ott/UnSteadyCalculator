function [p1,p2,p1up,p1low,p2up,p2low] = calc_isoline(SAMS,DEM,t,scenario,zm)
% Compute 10Be-14C isolines for step or spike erosion scenarios.
%   INPUTS
%   SAMS      Sample struct array (from cosmosampleread) with observed
%             concentrations and uncertainties (N10/N14, N10sigma/N14sigma).
%   DEM       GRIDobj used by cosmowatersheds to derive watershed properties
%             (lat/lon/elevation) for scaling factors and production rates.
%   t         Row/column vector of event ages (yr) where isolines are solved.
%   scenario  'step'  -> solve for [E1 E2] where E changes from E1 to E2.
%             'spike' -> solve for [E loss] where E is background erosion and
%                        loss is one-time stripped thickness.
%   zm        Soil-mixing depth in cm. Use zm = 0 for no mixing.
%
%   OUTPUTS
%   p1,p2     Best-fit parameters at each sample x time location.
%   p1up,p2up Best-fit parameters for (obs - sigma), used as one envelope.
%   p1low,p2low Best-fit parameters for (obs + sigma), used as one envelope.
%
%   NOTES
%   - Isoline inversion is non-convex. To avoid jumps to poor local minima
%     (for example E2/E1 spuriously collapsing toward zero), each time point
%     is solved with multiple starts and the lowest-misfit solution is kept.
%   - A warm-start from the previous time step is included to improve branch
%     continuity along each isoline.
%   - In 'spike' mode, curves are truncated at old ages once fitted loss
%     approaches the configured upper-loss ceiling. At that point most
%     inherited nuclides are effectively removed, so pre-spike history is
%     becomes unstable/non-unique.

% Relative mismatch threshold for spike isolines. 
misfitCutFrac = 0.02;  % percentage mismatch above which we assume solutions is useless

if nargin < 5 || isempty(zm)
    zm = 0;
end

SAMS = cosmowatersheds(SAMS,DEM);  % find catchments

lat = arrayfun(@(x) median(x.WSLat), SAMS);
lon = arrayfun(@(x) median(x.WSLon), SAMS);
alt = arrayfun(@(x) median(x.WSDEM.Z(:),'omitnan'), SAMS);

% catchment cosmo parameters
consts = make_constants();
spAll = sample_parameters(lat,lon,alt,consts);

% initiate arrays for calculation
n = numel(SAMS);
N10 = vertcat(SAMS.N10);
N14 = vertcat(SAMS.N14);
dN10 = vertcat(SAMS.N10sigma);
dN14 = vertcat(SAMS.N14sigma);

p1 = nan(n,length(t));
p2 = nan(n,length(t));
p1up = nan(n,length(t));
p2up = nan(n,length(t));
p1low = nan(n,length(t));
p2low = nan(n,length(t));

% set bounds for parameter search
switch scenario
    case 'step'
        x0 = [20, 300];
        UB = [1e3, 1e5];
    case 'spike'
        x0 = [50, 50];
        UB = [1e3, 1e4];
end
LB = [0, 0];

options = optimset('Display','off','MaxIter',2e4,'MaxFunEvals',2e4,'TolX',1e-6,'TolFun',1e-3);

for i = 1:n
    % Solve each catchment/sample independently.
    NlogicalRow = [~isnan(N10(i)), ~isnan(N14(i)), false];
    if sum(NlogicalRow(1:2)) < 2
        continue
    end

    obs = [N10(i); N14(i)];
    sig = [dN10(i); dN14(i)];

    % Extract just the current sample-specific parameters from the full set.
    sp = struct();
    fields = fieldnames(spAll);
    for k = 1:numel(fields)
        v = spAll.(fields{k});
        if isnumeric(v) && isvector(v) && numel(v) == numel(spAll.P10spal)
            sp.(fields{k}) = v(i);
        else
            sp.(fields{k}) = v;
        end
    end

    % Warm starts: keep continuity in t by seeding with previous solution.
    xPrev = x0;
    xPrevLow = x0;
    xPrevUp = x0;
    misfitBest = nan(1,length(t));
    misfitLow = nan(1,length(t));
    misfitUp = nan(1,length(t));

    for j = 1:length(t)
        tStep = [t(j), 0];

        % Combine default, warm, and a few exploratory seeds; keep best fit.
        switch scenario
            case 'step'
                starts = [x0; xPrev; 50, 5e3; 80, 3e4; 200, 10; 800, 1];
                startsLow = [x0; xPrevLow; 50, 5e3; 80, 3e4; 200, 10; 800, 1];
                startsUp = [x0; xPrevUp; 50, 5e3; 80, 3e4; 200, 10; 800, 1];
            case 'spike'
                starts = [x0; xPrev; 50, 500; 200, 100; 800, 20];
                startsLow = [x0; xPrevLow; 50, 500; 200, 100; 800, 20];
                startsUp = [x0; xPrevUp; 50, 500; 200, 100; 800, 20];
        end
        starts = unique(min(max(starts, LB), UB), 'rows', 'stable');
        startsLow = unique(min(max(startsLow, LB), UB), 'rows', 'stable');
        startsUp = unique(min(max(startsUp, LB), UB), 'rows', 'stable');

        % Best-fit isoline
        fun = @(x) misfit_model(x, scenario, tStep, sp, consts, zm, NlogicalRow, obs);
        [sol,exitflag] = optimize_multistart(fun, starts, LB, UB, options);
        if exitflag <= 0
            sol = [nan nan];
        else
            xPrev = sol;
        end
        p1(i,j) = sol(1); p2(i,j) = sol(2);
        if all(isfinite(sol))
            [~, predBest] = misfit_model(sol, scenario, tStep, sp, consts, zm, NlogicalRow, obs);
            misfitBest(j) = max(abs(predBest(:) - obs(:)) ./ max(abs(obs(:)), eps));
        end

        % Lower envelope (obs + sigma)
        targetLow = obs + sig;
        funLow = @(x) misfit_model(x, scenario, tStep, sp, consts, zm, NlogicalRow, targetLow);
        [solLow,exitflag] = optimize_multistart(funLow, startsLow, LB, UB, options);
        if exitflag <= 0
            solLow = [nan nan];
        else
            xPrevLow = solLow;
        end
        p1low(i,j) = solLow(1); p2low(i,j) = solLow(2);
        if all(isfinite(solLow))
            [~, predLow] = misfit_model(solLow, scenario, tStep, sp, consts, zm, NlogicalRow, targetLow);
            misfitLow(j) = max(abs(predLow(:) - targetLow(:)) ./ max(abs(targetLow(:)), eps));
        end

        % Upper envelope (obs - sigma)
        targetUp = obs - sig;
        funUp = @(x) misfit_model(x, scenario, tStep, sp, consts, zm, NlogicalRow, targetUp);
        [solUp,exitflag] = optimize_multistart(funUp, startsUp, LB, UB, options);
        if exitflag <= 0
            solUp = [nan nan];
        else
            xPrevUp = solUp;
        end
        p1up(i,j) = solUp(1); p2up(i,j) = solUp(2);
        if all(isfinite(solUp))
            [~, predUp] = misfit_model(solUp, scenario, tStep, sp, consts, zm, NlogicalRow, targetUp);
            misfitUp(j) = max(abs(predUp(:) - targetUp(:)) ./ max(abs(targetUp(:)), eps));
        end
    end

    if strcmp(scenario, 'spike')
        maxMisfit = max([misfitBest; misfitLow; misfitUp], [], 1);
        cutIdx = find(maxMisfit > misfitCutFrac, 1, 'first');
        if ~isempty(cutIdx)
            p1(i,cutIdx:end) = nan;
            p2(i,cutIdx:end) = nan;
            p1low(i,cutIdx:end) = nan;
            p2low(i,cutIdx:end) = nan;
            p1up(i,cutIdx:end) = nan;
            p2up(i,cutIdx:end) = nan;
        end
    end
end

% Clean obvious invalid solutions.
invalid = @(x) (~isfinite(x) | x < 0);
p1(invalid(p1)) = nan;   p2(invalid(p2)) = nan;
p1up(invalid(p1up)) = nan; p2up(invalid(p2up)) = nan;
p1low(invalid(p1low)) = nan; p2low(invalid(p2low)) = nan;

end


function [err, pred] = misfit_model(x, scenario, tStep, sp, consts, zm, NlogicalRow, obs)
switch scenario
    case 'step'
        pred = Nforward_discretized([x(1) x(2)], tStep, sp, consts, zm, 'step', NlogicalRow);
    case 'spike'
        pred = Nforward_discretized(x(1), tStep, sp, consts, zm, 'spike', NlogicalRow, x(2));
end
% L1 objective keeps the inversion robust to single-nuclide outliers.
err = sum(abs(pred(:) - obs(:)));
end


function [xBest,exitflagBest] = optimize_multistart(fun, starts, lb, ub, options)
xBest = [nan nan];
fBest = inf;
exitflagBest = -1;

for s = 1:size(starts,1)
    u0 = bounded_to_unbounded(starts(s,:), lb, ub);
    obj = @(u) fun(unbounded_to_bounded(u, lb, ub));
    [u,~,exitflag] = fminsearch(obj, u0, options);
    x = unbounded_to_bounded(u, lb, ub);

    if exitflag <= 0 || any(~isfinite(x))
        continue
    end

    fval = fun(x);
    if isfinite(fval) && fval < fBest
        fBest = fval;
        xBest = x;
        exitflagBest = exitflag;
    end
end
end


function u = bounded_to_unbounded(x, lb, ub)
z = (x - lb) ./ (ub - lb);
z = min(max(z,1e-10),1-1e-10);
u = log(z ./ (1 - z));
end


function x = unbounded_to_bounded(u, lb, ub)
s = 1 ./ (1 + exp(-u));
x = lb + (ub - lb) .* s;
end


