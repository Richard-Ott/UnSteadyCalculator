function [E, chg, Eup, Elow, chgup, chglow] = calc_curve_solver(SAMS, DEM, timebreaks, curvechanges)
% Solve for per-sample [background erosion, curve scale factor] from 10Be-14C.
%   INPUTS
%   SAMS          Sample struct array (from cosmosampleread) with observed
%                 concentrations and uncertainties (N10/N14, N10sigma/N14sigma).
%   DEM           GRIDobj used by cosmowatersheds to derive watershed properties.
%   timebreaks    Row vector of time points (yr BP, oldest to youngest, ending in 0).
%                 e.g. [10000, 6200, 700, 0]. The oldest interval sets the
%                 steady-state background; curve changes apply from timebreaks(2)
%                 onward. T passed to the forward model = timebreaks(2:end).
%   curvechanges  Relative erosion change per period relative to background.
%                 Length must equal length(timebreaks)-2.
%
%   OUTPUTS
%   E, chg        Best-fit background erosion (mm/ka) and curve scale factor, n x 1.
%   Eup, chgup    Parameter envelope from (obs - sigma), upper concentration bound.
%   Elow, chglow  Parameter envelope from (obs + sigma), lower concentration bound.

% calculate catchment sample parameters
SAMS = cosmowatersheds(SAMS, DEM);

lat = arrayfun(@(x) median(x.WSLat), SAMS);
lon = arrayfun(@(x) median(x.WSLon), SAMS);
alt = arrayfun(@(x) median(x.WSDEM.Z(:), 'omitnan'), SAMS);

consts = make_constants();
spAll = sample_parameters(lat, lon, alt, consts);

%% calculae erosion parameters

n = numel(SAMS);
N10  = vertcat(SAMS.N10);
N14  = vertcat(SAMS.N14);
dN10 = vertcat(SAMS.N10sigma);
dN14 = vertcat(SAMS.N14sigma);

E    = nan(n,1);  chg    = nan(n,1);
Eup  = nan(n,1);  chgup  = nan(n,1);
Elow = nan(n,1);  chglow = nan(n,1);

T  = timebreaks(2:end)';   % column vector — Nforward_discretized expects T' to be a row
% search boundaries
x0 = [50, 10];
LB = [0,   1];             % first parameter erosion, second change factor
UB = [1e3, 100];

% start at different points to make sure solver doesnt get stuck in local
% minima
starts = [x0; 20, 2; 200, 50; 500, 100; 50, 0.5; 100, 200];
starts = unique(min(max(starts, LB), UB), 'rows', 'stable');

options = optimset('Display','off','MaxIter',2e4,'MaxFunEvals',2e4,'TolX',1e-6,'TolFun',1e-3);

for i = 1:n
    NlogicalRow = [~isnan(N10(i)), ~isnan(N14(i)), false];
    if sum(NlogicalRow(1:2)) < 2
        continue
    end

    obs = [N10(i); N14(i)];
    sig = [dN10(i); dN14(i)];

    % Extract single-sample sp and attach the shared curvechange vector.
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
    sp.curvechange = curvechanges;

    % Best-fit solution
    fun = @(x) misfit_model(x, T, sp, consts, NlogicalRow, obs);
    [sol, exitflag] = optimize_multistart(fun, starts, LB, UB, options);
    if exitflag > 0
        E(i) = sol(1);  chg(i) = sol(2);
    end

    % Lower envelope (obs + sigma)
    funLow = @(x) misfit_model(x, T, sp, consts, NlogicalRow, obs + sig);
    [solLow, exitflag] = optimize_multistart(funLow, starts, LB, UB, options);
    if exitflag > 0
        Elow(i) = solLow(1);  chglow(i) = solLow(2);
    end

    % Upper envelope (obs - sigma)
    funUp = @(x) misfit_model(x, T, sp, consts, NlogicalRow, obs - sig);
    [solUp, exitflag] = optimize_multistart(funUp, starts, LB, UB, options);
    if exitflag > 0
        Eup(i) = solUp(1);  chgup(i) = solUp(2);
    end
end

end


function err = misfit_model(x, T, sp, consts, NlogicalRow, obs)
pred = Nforward_discretized(x(1), T, sp, consts, 'curve', NlogicalRow, x(2));
err  = sum(abs(pred(:) - obs(:)));
end


function [xBest, exitflagBest] = optimize_multistart(fun, starts, lb, ub, options)
xBest = [nan nan];
fBest = inf;
exitflagBest = -1;

for s = 1:size(starts, 1)
    u0 = bounded_to_unbounded(starts(s,:), lb, ub);
    obj = @(u) fun(unbounded_to_bounded(u, lb, ub));
    [u, ~, exitflag] = fminsearch(obj, u0, options);
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
z = min(max(z, 1e-10), 1-1e-10);
u = log(z ./ (1 - z));
end


function x = unbounded_to_bounded(u, lb, ub)
s = 1 ./ (1 + exp(-u));
x = lb + (ub - lb) .* s;
end
