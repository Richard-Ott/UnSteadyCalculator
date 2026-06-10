% This script computes the background erosion rate and curve scale factor
% for every sample using the pollen-based erosion curve from
% Jouffroy-Babiquot et al. 2021. Each sample is solved independently.
% Richard Ott, 2026

clc
clear
close all
addpath(genpath(pwd))

%% USER INPUT ----------------------------------------------------------- %

dataFile = 'data\WCdata_RFO.xlsx';
demFile  = '.\data\crete_clipped_utm.tif';

%% LOAD DATA ------------------------------------------------------------ %

SAMS = cosmosampleread(dataFile);
DEM  = GRIDobj(demFile);

%% BUILD CURVE ---------------------------------------------------------- %

curvedata  = load('./data/pollen.mat');
pollen     = curvedata.pollen;

timebreaks = [10000, 6200, 700, 0];   % yr BP; oldest interval = background
for i = 1:length(timebreaks) - 1
    timeRange    = pollen.yearsBP >= timebreaks(i+1) & pollen.yearsBP < timebreaks(i);
    meanPercTree(i) = mean(pollen.percTree(timeRange));
end
noTreePerc   = 100 - meanPercTree;
curvechanges = noTreePerc(2:end) ./ noTreePerc(1) - 1;   % relative changes vs background period

%% SOLVE PER SAMPLE ----------------------------------------------------- %

[E, chg, Eup, Elow, chgup, chglow] = calc_curve_solver(SAMS, DEM, timebreaks, curvechanges);

%% EROSION HISTORIES ---------------------------------------------------- %

erosion_histories = [E, E+ E.*chg.*curvechanges(1), E+ E.*chg.*curvechanges(2)]; 

%% PLOT ----------------------------------------------------------------- %

labels = {SAMS.ID};
nSamp  = numel(SAMS);
x      = 1:nSamp;

h1 = figure();
subplot(2,1,1)
errorbar(x, E, max(0, E - Elow), max(0, Eup - E), 'ko', 'LineWidth', 1.5, 'CapSize', 6)
set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45)
ylabel('Background erosion rate (mm/ka)')
title('Background erosion rate per sample')
grid on

subplot(2,1,2)
errorbar(x, chg, max(0, chg - chglow), max(0, chgup - chg), 'ko', 'LineWidth', 1.5, 'CapSize', 6)
set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45)
ylabel('Curve scale factor')
title('Curve scale factor per sample')
grid on

%% Pollen plot

% Sort so the fill patch is well behaved
[yearsSorted, idx] = sort(pollen.yearsBP(:), 'ascend');
percTreeSorted = pollen.percTree(idx);

valid = ~isnan(yearsSorted) & ~isnan(percTreeSorted);
yearsSorted = yearsSorted(valid);
percTreeSorted = percTreeSorted(valid);

% Dark, not-too-bright green
greenCol = [0.08 0.35 0.20];

h2 = figure('Color','w');
t = tiledlayout(2,1,'TileSpacing','tight','Padding','compact');

ax1 = nexttile;
hold(ax1,'on')

% Shade under pollen curve
fill([yearsSorted; flipud(yearsSorted)], ...
     [percTreeSorted; zeros(size(percTreeSorted))], ...
     greenCol, ...
     'FaceAlpha', 0.25, ...
     'EdgeColor', 'none');

% Main pollen line
plot(yearsSorted, percTreeSorted, ...
     'Color', greenCol, ...
     'LineWidth', 1.8);

% Mean tree cover as black step line over the time intervals
% stairs: y(i) is used from x(i) to x(i+1)
xStepPollen = [0 700 6200 10000];
yStepPollen = [meanPercTree(3) meanPercTree(2) meanPercTree(1) meanPercTree(1)];
stairs(xStepPollen, yStepPollen, ...
       'k', 'LineWidth', 2.0);

set(ax1, 'XDir', 'reverse', 'FontSize', 11, 'LineWidth', 1);
xlim(ax1, [0 10000]);
ylim(ax1, [0 max(100, ceil(max(percTreeSorted)/10)*10)]);
ylabel(ax1, 'Tree pollen (%)');
title(ax1, 'Pollen tree cover through time');
box(ax1,'on')

% Bottom subplot: erosion
ax2 = nexttile;
hold(ax2,'on')

% Model erosion history:
xStepErosion = [0 700 6300 10000];
yStepErosion = [mean(erosion_histories(:,3)) mean(erosion_histories(:,2))...
    mean(erosion_histories(:,1)) mean(erosion_histories(:,1))];

hold on
for i =1:size(erosion_histories,1)
    stairs(xStepErosion, [fliplr(erosion_histories(i,:)),erosion_histories(i,1) ], ...
      'Color', [.5, .5, .5], 'LineWidth', 1.0);
end

stairs(xStepErosion, yStepErosion, ...
       'k', 'LineWidth', 2.2);

set(ax2, 'XDir', 'reverse', 'FontSize', 11, 'LineWidth', 1);
xlim(ax2, [0 10000]);
% ylim(ax2, [0 1000]);
xlabel(ax2, 'Years BP');
ylabel(ax2, 'Erosion rate (mm kyr^{-1})');
title(ax2, 'Modelled erosion through time');
box(ax2,'on')
grid(ax2,'on')

% Keep x-axes aligned
linkaxes([ax1 ax2], 'x')
%% export

set(h2, 'PaperOrientation', 'portrait')
print(h2, 'Curve_solver_WC.pdf', '-dpdf', '-vector', '-bestfit')