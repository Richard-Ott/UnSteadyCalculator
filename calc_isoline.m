function [p1,p2,p1up,p1low,p2up,p2low] = calc_isoline(SAMS,DEM,t,scenario)
%CALC_ISOLINE Compute 10Be-14C isolines for step or spike erosion scenarios.
%
%   [p1,p2,p1up,p1low,p2up,p2low] = calc_isoline(SAMS,DEM,t,scenario)

if nargin < 4
    error('calc_isoline requires SAMS, DEM, t, and scenario.');
end
if ~ismember(scenario, {'step','spike'})
    error('scenario must be ''step'' or ''spike''.');
end

SAMS = cosmowatersheds(SAMS,DEM);

lat = arrayfun(@(x) median(x.WSLat), SAMS);
lon = arrayfun(@(x) median(x.WSLon), SAMS);
alt = arrayfun(@(x) median(x.WSDEM.Z(:),'omitnan'), SAMS);

consts = make_constants();
spAll = sample_parameters(lat,lon,alt,consts);

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
spikeLossCutFrac = 0.90;

for i = 1:n
    NlogicalRow = [~isnan(N10(i)), ~isnan(N14(i)), false];
    if sum(NlogicalRow(1:2)) < 2
        continue
    end

    obs = [N10(i); N14(i)];
    sig = [dN10(i); dN14(i)];

    sp = slice_sample_parameters(spAll, i);
    xPrev = x0;
    xPrevLow = x0;
    xPrevUp = x0;

    for j = 1:length(t)
        tStep = [t(j), 0];

        fun = @(x) misfit_model(x, scenario, tStep, sp, consts, NlogicalRow, obs);
        starts = build_starts(scenario, x0, xPrev, LB, UB);
        [sol,exitflag] = optimize_multistart(fun, starts, LB, UB, options);
        if exitflag <= 0
            sol = [nan nan];
        else
            xPrev = sol;
        end
        p1(i,j) = sol(1); p2(i,j) = sol(2);

        funLow = @(x) misfit_model(x, scenario, tStep, sp, consts, NlogicalRow, obs + sig);
        startsLow = build_starts(scenario, x0, xPrevLow, LB, UB);
        [solLow,exitflag] = optimize_multistart(funLow, startsLow, LB, UB, options);
        if exitflag <= 0
            solLow = [nan nan];
        else
            xPrevLow = solLow;
        end
        p1low(i,j) = solLow(1); p2low(i,j) = solLow(2);

        funUp = @(x) misfit_model(x, scenario, tStep, sp, consts, NlogicalRow, obs - sig);
        startsUp = build_starts(scenario, x0, xPrevUp, LB, UB);
        [solUp,exitflag] = optimize_multistart(funUp, startsUp, LB, UB, options);
        if exitflag <= 0
            solUp = [nan nan];
        else
            xPrevUp = solUp;
        end
        p1up(i,j) = solUp(1); p2up(i,j) = solUp(2);
    end

    if strcmp(scenario, 'spike')
        cutIdx = detect_spike_cutoff(p2(i,:), p2low(i,:), p2up(i,:), UB(2), spikeLossCutFrac);
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

invalid = @(x) (~isfinite(x) | x < 0);
p1(invalid(p1)) = nan;   p2(invalid(p2)) = nan;
p1up(invalid(p1up)) = nan; p2up(invalid(p2up)) = nan;
p1low(invalid(p1low)) = nan; p2low(invalid(p2low)) = nan;

end


function err = misfit_model(x, scenario, tStep, sp, consts, NlogicalRow, obs)
switch scenario
    case 'step'
        pred = Nforward_discretized([x(1) x(2)], tStep, sp, consts, 'step', NlogicalRow);
    case 'spike'
        pred = Nforward_discretized(x(1), tStep, sp, consts, 'spike', NlogicalRow, x(2));
end
err = sum(abs(pred(:) - obs(:)));
end


function [x,exitflag] = optimize_bounded(fun, x0, lb, ub, options)
u0 = bounded_to_unbounded(x0, lb, ub);
obj = @(u) fun(unbounded_to_bounded(u, lb, ub));
[u,~,exitflag] = fminsearch(obj, u0, options);
x = unbounded_to_bounded(u, lb, ub);
end


function [xBest,exitflagBest] = optimize_multistart(fun, starts, lb, ub, options)
xBest = [nan nan];
fBest = inf;
exitflagBest = -1;

for s = 1:size(starts,1)
    x0 = starts(s,:);
    [x,exitflag] = optimize_bounded(fun, x0, lb, ub, options);
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


function starts = build_starts(scenario, xDefault, xPrev, lb, ub)
switch scenario
    case 'step'
        base = [
            xDefault
            xPrev
            50, 5e3
            80, 3e4
            200, 10
            800, 1
        ];
    case 'spike'
        base = [
            xDefault
            xPrev
            50, 500
            200, 100
            800, 20
        ];
end

starts = max(base, repmat(lb, size(base,1), 1));
starts = min(starts, repmat(ub, size(base,1), 1));
starts = unique(starts, 'rows', 'stable');
end


function cutIdx = detect_spike_cutoff(lossBest, lossLow, lossUp, lossUB, cutFrac)
lossThresh = cutFrac * lossUB;
lossHighForUnstable = 0.50 * lossUB;
jumpDexThresh = 0.35;

highLoss = (lossBest >= lossThresh) | (lossLow >= lossThresh) | (lossUp >= lossThresh);
highLoss = highLoss & isfinite(lossBest);

unstableHighLoss = false(size(lossBest));
for j = 3:numel(lossBest)
    if ~all(isfinite(lossBest(j-2:j)))
        continue
    end
    l0 = max(lossBest(j-2), 1e-8);
    l1 = max(lossBest(j-1), 1e-8);
    l2 = max(lossBest(j), 1e-8);
    d1 = abs(log10(l1) - log10(l0));
    d2 = abs(log10(l2) - log10(l1));
    zigzag = sign(lossBest(j) - lossBest(j-1)) ~= sign(lossBest(j-1) - lossBest(j-2));
    unstableHighLoss(j) = zigzag && ((d1 > jumpDexThresh) || (d2 > jumpDexThresh)) && (lossBest(j) >= lossHighForUnstable);
end

trigger = highLoss | unstableHighLoss;

cutIdx = [];
if isempty(trigger)
    return
end

for j = 1:(numel(trigger)-1)
    if trigger(j) && trigger(j+1)
        cutIdx = j;
        return
    end
end

if trigger(end)
    cutIdx = numel(trigger);
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


function sp = slice_sample_parameters(spAll, idx)
sp = struct();
f = fieldnames(spAll);
for k = 1:numel(f)
    v = spAll.(f{k});
    if isnumeric(v) && (isvector(v) && numel(v) == numel(spAll.P10spal))
        sp.(f{k}) = v(idx);
    else
        sp.(f{k}) = v;
    end
end
end