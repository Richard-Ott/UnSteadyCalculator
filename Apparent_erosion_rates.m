% This script calculates apparent erosion rates for 10Be, 14C, and 26Al
% with the Nforward_steadystate function.
% Richard Ott, 2026

clc
clear
close all
addpath(genpath(pwd))

%% USER INPUT ----------------------------------------------------------- %

demFile = '.\data\crete_clipped_utm.tif';
dataFile = 'data\WCdata_RFO.xlsx';

emin = 1;      % minimum erosion rate searched (mm/ka)
emax = 1e4;       % maximum erosion rate searched (mm/ka)

%% LOAD DATA ------------------------------------------------------------ %

DEM = GRIDobj(demFile);
SAMS = cosmosampleread(dataFile);
SAMS = cosmowatersheds(SAMS, DEM);

data = readtable(dataFile);
consts = make_constants();

lat = arrayfun(@(x) median(x.WSLat), SAMS);
lon = arrayfun(@(x) median(x.WSLon), SAMS);
alt = arrayfun(@(x) median(x.WSDEM.Z(:), 'omitnan'), SAMS);
sp = sample_parameters(lat, lon, alt, consts);

n = numel(SAMS);
samp_names = {SAMS.ID}';

Nlogical = [~isnan(data.N10), ~isnan(data.N14), ~isnan(data.N26)];  
Nobs = [data.N10, data.N14, data.N26];
dNobs = [data.N10sigma, data.N14sigma, data.N26sigma];

%% CALCULATE APPARENT STEADY-STATE EROSION RATES ----------------------- %

Erate10 = nan(n,1);
dErate10 = nan(n,1);
Erate14 = nan(n,1);
dErate14 = nan(n,1);
Erate26 = nan(n,1);
dErate26 = nan(n,1);

%% Calculate 10Be E-rates
for i = 1:n
    if Nlogical(i,1)
        logE = fminbnd(@(x) steady_state_misfit(x, Nobs(i,1), i, 1, sp, consts), log10(emin), log10(emax));
        Erate10(i) = 10.^logE;

        logEhigh = fminbnd(@(x) steady_state_misfit(x, Nobs(i,1) - dNobs(i,1), i, 1, sp, consts), log10(emin), log10(emax));
        logElow = fminbnd(@(x) steady_state_misfit(x, Nobs(i,1) + dNobs(i,1), i, 1, sp, consts), log10(emin), log10(emax));
        dErate10(i) = 0.5 * abs(10.^logEhigh - 10.^logElow);
    end
end

%% Calculate 14C E-rates
for i = 1:n
    if Nlogical(i,2)
        logE = fminbnd(@(x) steady_state_misfit(x, Nobs(i,2), i, 2, sp, consts), log10(emin), log10(emax));
        Erate14(i) = 10.^logE;

        logEhigh = fminbnd(@(x) steady_state_misfit(x, Nobs(i,2) - dNobs(i,2), i, 2, sp, consts), log10(emin), log10(emax));
        logElow = fminbnd(@(x) steady_state_misfit(x, Nobs(i,2) + dNobs(i,2), i, 2, sp, consts), log10(emin), log10(emax));
        dErate14(i) = 0.5 * abs(10.^logEhigh - 10.^logElow);
    end
end

%% Calculate 26Al E-rates
for i = 1:n
    if Nlogical(i,3)
        logE = fminbnd(@(x) steady_state_misfit(x, Nobs(i,3), i, 3, sp, consts), log10(emin), log10(emax));
        Erate26(i) = 10.^logE;

        logEhigh = fminbnd(@(x) steady_state_misfit(x, Nobs(i,3) - dNobs(i,3), i, 3, sp, consts), log10(emin), log10(emax));
        logElow = fminbnd(@(x) steady_state_misfit(x, Nobs(i,3) + dNobs(i,3), i, 3, sp, consts), log10(emin), log10(emax));
        dErate26(i) = 0.5 * abs(10.^logEhigh - 10.^logElow);
    end
end

%% OUTPUT --------------------------------------------------------------- %

vars = {'Sample', 'E10_mmka', 'dE10_mmka', 'E14_mmka', 'dE14_mmka', 'E26_mmka', 'dE26_mmka'};
out_table = table(samp_names, Erate10, dErate10, Erate14, dErate14, Erate26, dErate26, 'VariableNames', vars);
disp(out_table)

%% PLOT ----------------------------------------------------------------- %
labelOffset = 10; % label offset in plot units

C1 = [0, 92/255, 171/255];
short_names = cellfun(@(s) regexprep(regexp(s, '(?<=-)[^-]+$', 'match', 'once'), '''$', ''), ...
    samp_names, 'UniformOutput', false);

has14 = any(isfinite(Erate10) & isfinite(Erate14));
has26 = any(isfinite(Erate10) & isfinite(Erate26));

if has14 || has26
    nPanels = has14 + has26;
    figure()
    panelIdx = 1;

    if has14
        valid = isfinite(Erate10) & isfinite(Erate14);
        subplot(1, nPanels, panelIdx)
        errorbar(Erate10(valid), Erate14(valid), dErate14(valid), dErate14(valid), dErate10(valid), dErate10(valid), 'o', ...
            'MarkerSize', 8, 'MarkerEdgeColor', [.3, .3, .3], 'MarkerFaceColor', C1, ...
            'CapSize', 0, 'Color', [.3, .3, .3])
        hold on
        text(Erate10(valid) + labelOffset / 3, Erate14(valid) + labelOffset * 4, short_names(valid))
        plot(0:1:2000, 0:1:2000, '--k')
        plot(0:1:2000, 0:5:10000, '--k')
        plot(0:1:2000, 0:10:20000, '--k')
        xlim([0, 320])
        ylim([0, 1500])
        axis square
        xlabel('$\epsilon_{app}\,^{10}\mathrm{Be}$ mm/ka', 'Interpreter', 'latex')
        ylabel('$\epsilon_{app}\,^{14}\mathrm{C}$ mm/ka', 'Interpreter', 'latex')
        title('^{14}C vs ^{10}Be', 'Interpreter', 'tex')
        panelIdx = panelIdx + 1;


        % 1-1 1-5 and 1-10 text label
        % get axis limits
        xl = xlim;
        yl = ylim;
        
        % compute angle correction so text aligns with lines on screen
        pbaspect([1 1 1]) % because you use axis square
        angle1 = atan2d(1*(xl(2)-xl(1)),(yl(2)-yl(1)));
        angle5 = atan2d(5*(xl(2)-xl(1)),(yl(2)-yl(1)));
        angle10 = atan2d(10*(xl(2)-xl(1)),(yl(2)-yl(1)));
        
        % label positions (chosen inside the visible axes)
        x1 = 300;  y1 = 300;
        x5 = 270;   y5 = 5*x5;
        x10 =130;  y10 = 10*x10;
        
        text(x1,y1,'1:1','Rotation',angle1,'HorizontalAlignment','left','FontWeight','bold')
        text(x5,y5,'5:1','Rotation',angle5,'HorizontalAlignment','left','FontWeight','bold')
        text(x10,y10,'10:1','Rotation',angle10,'HorizontalAlignment','left','FontWeight','bold')
    end

    if has26
        valid = isfinite(Erate10) & isfinite(Erate26);
        subplot(1, nPanels, panelIdx)
        errorbar(Erate10(valid), Erate26(valid), dErate26(valid), dErate26(valid), dErate10(valid), dErate10(valid), 'o', ...
            'MarkerSize', 8, 'MarkerEdgeColor', [.3, .3, .3], 'MarkerFaceColor', C1, ...
            'CapSize', 0, 'Color', [.3, .3, .3])
        hold on
        text(Erate10(valid) + labelOffset / 3, Erate26(valid) + labelOffset * 4, short_names(valid))
        plot(0:1:2000, 0:1:2000, '--k')
        plot(0:1:2000, 0:5:10000, '--k')
        plot(0:1:2000, 0:10:20000, '--k')
        xlim([0, 320])
        ylim([0, 1500])
        axis square
        xlabel('$\epsilon_{app}\,^{10}\mathrm{Be}$ mm/ka', 'Interpreter', 'latex')
        ylabel('$\epsilon_{app}\,^{26}\mathrm{Al}$ mm/ka', 'Interpreter', 'latex')
        title('^{26}Al vs ^{10}Be', 'Interpreter', 'tex')
    end
end

% forward model call
function err = steady_state_misfit(logE, Ntarget, sampleIdx, nuclideIdx, sp, consts)
Nlogical = false(numel(sp.P10spal), 3);
Nlogical(sampleIdx, nuclideIdx) = true;
Nss = Nforward_steadystate(10.^logE, sp, consts, Nlogical);
err = (log(Nss) - log(Ntarget)).^2;
end
