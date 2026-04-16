function [E10,E14,E10up,E10low,E14up,E14low,Ebest,zmbest,hasSolution,fitMisfit] = calc_mixing_depth(SAMS,DEM,zm)
%CALC_MIXING_DEPTH Search for steady-state erosion rate and mixing depth.
%
%   [E10,E14,E10up,E10low,E14up,E14low,Ebest,zmbest,hasSolution,fitMisfit] = ...
%       calc_mixing_depth(SAMS,DEM,zm)
%
%   This function  assumes constant
%   erosion through time and asks the following for each sample:
%
%       Which erosion rate E and which mixing depth z_m reproduce the
%       measured 10Be and 14C concentrations at the same time?
%
%   The workflow is:
%     1) For every candidate z_m in the input vector, solve the apparent
%        steady-state erosion rate from 10Be alone and from 14C alone.
%     2) Compare those two curves E10(z_m) and E14(z_m).
%     3) Use their closest approach as the start point for a joint search in
%        [E, z_m] space.
%     4) Keep the solution only if both nuclides are reproduced within a
%        small tolerance.
%
%   INPUTS
%   SAMS   Sample struct array from cosmosampleread.
%   DEM    GRIDobj used by cosmowatersheds.
%   zm     Vector of candidate mixing depths in cm.
%
%   OUTPUTS
%   E10,E14      Apparent erosion-rate curves (mm/ka) vs z_m for each sample.
%   E10up/low    10Be curves using N10 -/+ sigma.
%   E14up/low    14C curves using N14 -/+ sigma.
%   Ebest        Best-fit common erosion rate per sample (mm/ka).
%   zmbest       Best-fit common mixing depth per sample (cm).
%   hasSolution  Logical flag indicating whether a joint solution was found.
%   fitMisfit    Final sigma-weighted misfit of the joint solution.
%
% Richard Ott, 2026

% acceptance threshold for found solution
thres = 0.02; % 0.02 = 2% within observed nuclide concentration

if nargin < 3 || isempty(zm)
    zm = 0:10:300;
end
zm = zm(:)';   % work internally with a row vector for easier plotting

% Derive catchment geometry and all production parameters once. Those values
% are then reused for each sample-specific inversion below.
SAMS = cosmowatersheds(SAMS,DEM);
lat = arrayfun(@(x) median(x.WSLat), SAMS);
lon = arrayfun(@(x) median(x.WSLon), SAMS);
alt = arrayfun(@(x) median(x.WSDEM.Z(:),'omitnan'), SAMS);

consts = make_constants();
spAll = sample_parameters(lat,lon,alt,consts);

n = numel(SAMS);
m = numel(zm);

% Observations and 1-sigma analytical uncertainty for the two nuclides used
% in the mixing-depth search.
N10 = vertcat(SAMS.N10);
N14 = vertcat(SAMS.N14);
dN10 = vertcat(SAMS.N10sigma);
dN14 = vertcat(SAMS.N14sigma);

% Preallocate outputs. Each row corresponds to one sample and each column
% corresponds to one candidate mixing depth in zm.
E10 = nan(n,m);
E14 = nan(n,m);
E10up = nan(n,m);
E10low = nan(n,m);
E14up = nan(n,m);
E14low = nan(n,m);
Ebest = nan(n,1);
zmbest = nan(n,1);
hasSolution = false(n,1);
fitMisfit = nan(n,1);

% Search bounds for the steady-state solution.
eBounds = [1, 1e4];
zBounds = [min(zm), max(zm)];

% Separate settings for the 1-D and 2-D searches.
opt1 = optimset('Display','off','MaxIter',5e3,'MaxFunEvals',5e3,'TolX',1e-6,'TolFun',1e-8);
opt2 = optimset('Display','off','MaxIter',2e4,'MaxFunEvals',2e4,'TolX',1e-6,'TolFun',1e-8);

for i = 1:n
    % Extract the production parameters for one sample 
    sp = struct();
    f = fieldnames(spAll);
    for k = 1:numel(f)
        v = spAll.(f{k});
        if isnumeric(v) && isvector(v) && numel(v) == numel(spAll.P10spal)
            sp.(f{k}) = v(i);
        else
            sp.(f{k}) = v;
        end
    end

    % ------------------------------------------------------------------
    % 1) Build E(z_m) curves from each nuclide separately.
    % ------------------------------------------------------------------
    for j = 1:m
        zHere = zm(j);
        E10(i,j) = solve_single_nuclide(N10(i), 1, sp, consts, zHere, eBounds, opt1);
        E14(i,j) = solve_single_nuclide(N14(i), 2, sp, consts, zHere, eBounds, opt1);

        % Repeat the same search for obs ± sigma to provide an uncertainty
        % envelope around the apparent erosion-rate curves.
        E10up(i,j) = solve_single_nuclide(max(N10(i) - dN10(i), eps), 1, sp, consts, zHere, eBounds, opt1);
        E10low(i,j) = solve_single_nuclide(N10(i) + dN10(i), 1, sp, consts, zHere, eBounds, opt1);
        E14up(i,j) = solve_single_nuclide(max(N14(i) - dN14(i), eps), 2, sp, consts, zHere, eBounds, opt1);
        E14low(i,j) = solve_single_nuclide(N14(i) + dN14(i), 2, sp, consts, zHere, eBounds, opt1);
    end

    % ------------------------------------------------------------------
    % 2) Use the full 10Be and 14C curves directly.
    % ------------------------------------------------------------------
    dCurve = log10(E10(i,:)) - log10(E14(i,:));
    eValid = sqrt(E10(i,:) .* E14(i,:));
    [~, idx0] = min(abs(dCurve));

    xDefault = [eValid(idx0), zm(idx0)];
    starts = [
        xDefault
        20, 0
        50, min(25, zBounds(2))
        100, min(50, zBounds(2))
        200, min(100, zBounds(2))
        500, min(200, zBounds(2))
    ];
    starts(:,1) = min(max(starts(:,1), eBounds(1)), eBounds(2));
    starts(:,2) = min(max(starts(:,2), zBounds(1)), zBounds(2));
    starts = unique(starts, 'rows', 'stable');

    obs = [N10(i); N14(i)];
    sig = [dN10(i); dN14(i)];

    % ------------------------------------------------------------------
    % 3) Jointly optimize the common erosion rate and mixing depth.
    %    The objective is written in sigma-scaled log space so that both
    %    nuclides contribute comparably to the fit.
    % ------------------------------------------------------------------
    fun = @(x) concentration_misfit(x(1), x(2), obs, sig, [true true false], sp, consts);
    [xbest, exitflag] = optimize_multistart(fun, starts, [eBounds(1), zBounds(1)], [eBounds(2), zBounds(2)], opt2);

    if exitflag <= 0 || any(~isfinite(xbest))
        continue
    end

    pred = Nforward_steadystate(xbest(1), sp, consts, xbest(2), [true true false]);
    relErr = abs(pred(:) - obs) ./ max(obs, eps);
    fitMisfit(i) = fun(xbest);

    % ------------------------------------------------------------------
    % 4) Accept the solution only if the two measured concentrations can be
    %    matched reasonably well.
    % ------------------------------------------------------------------
    if all(relErr < thres) 
        Ebest(i) = xbest(1);
        zmbest(i) = xbest(2);
        hasSolution(i) = true;
    end
end

end


function Ebest = solve_single_nuclide(Ntarget, nuclideIdx, sp, consts, zm, eBounds, options)
% Solve for the steady-state erosion rate that matches one nuclide at a
% fixed mixing depth. The search is done in log(E) space because the valid
% erosion-rate range spans several orders of magnitude.

if ~isfinite(Ntarget) || Ntarget <= 0
    Ebest = nan;
    return
end

logLB = log10(eBounds(1));
logUB = log10(eBounds(2));
seedLogE = unique([mean([logLB logUB]), log10([5, 20, 100, 500, 2000])]);

bestF = inf;
Ebest = nan;
for k = 1:numel(seedLogE)
    x0 = min(max(seedLogE(k), logLB), logUB);
    u0 = bounded_to_unbounded(x0, logLB, logUB);

    Nlogical = false(1,3);
    Nlogical(nuclideIdx) = true;
    obj = @(u) concentration_misfit(10.^unbounded_to_bounded(u, logLB, logUB), zm, Ntarget, [], Nlogical, sp, consts);
    [u, fval, exitflag] = fminsearch(obj, u0, options);

    if exitflag > 0 && isfinite(fval) && (fval < bestF)
        bestF = fval;
        Ebest = 10.^unbounded_to_bounded(u, logLB, logUB);
    end
end
end


function [xBest,exitflagBest] = optimize_multistart(fun, starts, lb, ub, options)
% Run several fminsearch starts while enforcing bounds through a logistic
% transform. 

xBest = nan(1, size(starts,2));
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


function err = concentration_misfit(E, zm, obs, sig, Nlogical, sp, consts)
% Shared misfit helper used by both the 1-D single-nuclide search and the
% final 2-D joint [E, z_m] search.
pred = Nforward_steadystate(E, sp, consts, zm, Nlogical);

obs = obs(:);
pred = pred(:);
if any(~isfinite(pred)) || any(pred <= 0) || any(~isfinite(obs)) || any(obs <= 0)
    err = inf;
    return
end

if isempty(sig)
    % For the single-nuclide case, plain log-space distance is enough.
    err = sum(abs(log(pred) - log(obs)));
else
    % For the joint search, weight by analytical uncertainty so both
    % nuclides contribute in a balanced way.
    sigEff = max(sig(:), 0.01 * obs);
    relSigma = sigEff ./ max(obs, eps);
    err = sum(((log(pred) - log(obs)) ./ relSigma).^2);
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



