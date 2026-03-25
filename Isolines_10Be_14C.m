% This code computes the isolines of erosion ratios for a step or spike
% change in erosion rate, given 10Be-14C nuclice concentration measurements.
% The code might take some time to run, mainly because the Matlab
% optimization function can get stuck in some local minima. Thus, I run
% several optimizations and only plot the results from the best one.
% Especially for the spike scenario, at high spikes you might have removed
% all of the nuclides, and thus, there's nothing to know about the nuclide
% concentrations beyond that point. 
% Richard Ott, 2024

clc
clear
close all
addpath(genpath(pwd))

%% USER INPUT ----------------------------------------------------------- %

scenario = 'step'; % 'step' or 'spike'

dataFile = 'data\WCdata_RFO.xlsx';
demFile = '.\data\crete_clipped_utm.tif';

t = 1e2:2e2:1e4;              % years since step/spike event

%% LOAD DATA ------------------------------------------------------------ %

SAMS = cosmosampleread(dataFile);
DEM = GRIDobj(demFile);

%% CALCULATE ISOLINES --------------------------------------------------- %

[p1,p2,p1up,p1low,p2up,p2low] = calc_isoline(SAMS, DEM, t, scenario);

% for spike scenarios, one often reaches apoint where the spike is huge
% (>1500) and erases almost the entire production profile. From that point
% back in time, you won;t be able to constrain a solution. Here I
% discard all of these values
if strcmp(scenario, 'spike')
    maxLossCm = 1000;
    for i = 1:size(p2,1)
        exceed = (p2(i,:) > maxLossCm) | (p2up(i,:) > maxLossCm) | (p2low(i,:) > maxLossCm);
        firstBad = find(exceed, 1, 'first');
        if ~isempty(firstBad)
            p1(i,firstBad:end) = nan;
            p2(i,firstBad:end) = nan;
            p1up(i,firstBad:end) = nan;
            p2up(i,firstBad:end) = nan;
            p1low(i,firstBad:end) = nan;
            p2low(i,firstBad:end) = nan;
        end
    end
end

%% PLOT ----------------------------------------------------------------- %

cc = parula(numel(SAMS));

switch scenario
    case 'step'
        labels = {SAMS.ID};
        figure()
        for i = 1:numel(SAMS)
            ratio = p2(i,:) ./ p1(i,:);
            ratioUp = p2up(i,:) ./ p1up(i,:);
            ratioLow = p2low(i,:) ./ p1low(i,:);

            plot(t, ratio, '-', 'Color', cc(i,:), 'LineWidth', 1.5)
            hold on

            indsUp = isfinite(ratioUp);
            indsLow = isfinite(ratioLow);
            if any(indsUp) && any(indsLow)
                f = fill([t(indsUp), fliplr(t(indsLow))], ...
                    [ratioUp(indsUp), fliplr(ratioLow(indsLow))], ...
                    cc(i,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
                set(get(get(f,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
        end

        xlabel('years since step change')
        ylabel('E2 / E1')
        title('Step-change isolines from 10Be-14C')
        legend(labels, 'Interpreter', 'none', 'Location', 'best')

    case 'spike'
        figure()
        for i = 1:numel(SAMS)
            subplot(1,2,1)
            plot(t, p2(i,:), '-', 'Color', cc(i,:), 'LineWidth', 1.5)
            hold on

            subplot(1,2,2)
            plot(t, p1(i,:), '-', 'Color', cc(i,:), 'LineWidth', 1.5)
            hold on
        end

        subplot(1,2,1)
        xlabel('years since spike')
        ylabel('Loss (cm)')
        title('Spike isolines: loss')

        subplot(1,2,2)
        xlabel('years since spike')
        ylabel('Background erosion E (mm/ka)')
        title('Spike isolines: erosion')

        labels = {SAMS.ID};
        legend(labels, 'Interpreter', 'none', 'Location', 'best')
end

% ylim([0 350])
