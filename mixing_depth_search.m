% This script searches for the steady-state erosion rate and soil mixing
% depth that jointly explain the measured 10Be and 14C concentrations for
% each sample. It also plots erosion-rate versus mixing-depth curves and
% highlights samples where no common solution is found.
% Richard Ott, 2026

clc
clear
close all
addpath(genpath(pwd))

%% USER INPUT ----------------------------------------------------------- %

demFile = '.\data\crete_clipped_utm.tif';
dataFile = 'data\WCdata_RFO.xlsx';

zm = 0:5:500;   % search vector for soil mixing depth in cm

%% LOAD DATA ------------------------------------------------------------ %

SAMS = cosmosampleread(dataFile);
DEM = GRIDobj(demFile);

%% CALCULATE MIXING-DEPTH SOLUTIONS ------------------------------------ %

[E10,E14,E10up,E10low,E14up,E14low,Ebest,zmbest,hasSolution,fitMisfit] = calc_mixing_depth(SAMS, DEM, zm);

labels = {SAMS.ID}';
short_names = cellfun(@(s) regexprep(regexp(s, '(?<=-)[^-]+$', 'match', 'once'), '''$', ''), ...
    labels, 'UniformOutput', false);

out_table = table(labels, hasSolution, Ebest, zmbest, fitMisfit, ...
    'VariableNames', {'Sample', 'HasSolution', 'E_mmka', 'zm_cm', 'Misfit'});
disp(out_table)

%% PLOT CURVES FOR EACH SAMPLE ----------------------------------------- %

n = numel(SAMS);
nCols = min(3, max(1, ceil(sqrt(n))));
nRows = ceil(n / nCols);

C10 = [0, 92/255, 171/255];
C14 = [190/255, 30/255, 45/255];

figure('Color', 'w')
tl = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Steady-state erosion versus soil mixing depth from ^{10}Be and ^{14}C')

h10 = gobjects(1,1);
h14 = gobjects(1,1);
hFit = gobjects(1,1);

for i = 1:n
    nexttile
    hold on

    idx10 = isfinite(E10(i,:));
    idx14 = isfinite(E14(i,:));
    idx10band = isfinite(E10up(i,:)) & isfinite(E10low(i,:));
    idx14band = isfinite(E14up(i,:)) & isfinite(E14low(i,:));

    if any(idx10band)
        fill([zm(idx10band), fliplr(zm(idx10band))], ...
             [E10up(i,idx10band), fliplr(E10low(i,idx10band))], ...
             C10, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end
    if any(idx14band)
        fill([zm(idx14band), fliplr(zm(idx14band))], ...
             [E14up(i,idx14band), fliplr(E14low(i,idx14band))], ...
             C14, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    end

    if any(idx10)
        h10 = plot(zm, E10(i,:), '-', 'Color', C10, 'LineWidth', 1.5);
    end
    if any(idx14)
        h14 = plot(zm, E14(i,:), '-', 'Color', C14, 'LineWidth', 1.5);
    end

    if hasSolution(i)
        hFit = plot(zmbest(i), Ebest(i), 'ko', 'MarkerFaceColor', [1 .85 0], 'MarkerSize', 6);
        title(short_names{i}, 'Interpreter', 'none')
    else
        title([short_names{i} '  no solution'], 'Interpreter', 'none', 'Color', [0.75 0 0])
        text(0.04, 0.92, 'no common E-z_m solution', 'Units', 'normalized', ...
            'Color', [0.75 0 0], 'FontSize', 8, 'FontWeight', 'bold');
    end

    xlabel('z_m (cm)')
    ylabel('\epsilon (mm/ka)')
    xlim([min(zm), max(zm)])
    grid on
    box on
end

legend([h10 h14 hFit], {'^{10}Be', '^{14}C', 'joint solution'}, 'Location', 'bestoutside')

%% SUMMARY PLOT --------------------------------------------------------- %

figure('Color', 'w')
hold on
solved = hasSolution & isfinite(Ebest) & isfinite(zmbest);

if any(solved)
    scatter(zmbest(solved), Ebest(solved), 45, 'filled', 'MarkerFaceColor', 'k')
    text(zmbest(solved) + 2, Ebest(solved), short_names(solved), 'FontSize', 8)
end

grid on
box on
xlabel('z_m (cm)')
ylabel('\epsilon (mm/ka)')
title('Best-fit steady-state erosion rate and mixing depth per sample')

if any(~solved)
    unsolvedNames = strjoin(short_names(~solved), newline);
    annotation('textbox', [0.64 0.68 0.3 0.22], 'String', ['No solution found:' newline unsolvedNames], ...
        'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', [0.75 0 0], 'Color', [0.75 0 0]);
end
