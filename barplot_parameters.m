function H=barplot_parameters(m,varnames,prior_ranges,varargin)
% Plots grouped MCMC chain parameters using separate boxplots per variable group
% Additional arguments:
%   'bestmodel' - best model to overlay
%   'truevals'  - true parameter values for reference

% Richard Ott, 2025

% Parse optional inputs
p = inputParser;
addParameter(p, 'bestmodel', []);
addParameter(p, 'truevals', []);
parse(p, varargin{:});

% Gather all accepted samples
models = m(:,:); % collapse all model dimensions


group_prefixes = cell(size(varnames));
for i = 1:numel(varnames)
    group_prefixes{i} = regexp(varnames{i}, '^[^\d]*', 'match', 'once');
end



[unique_groups, ~, group_ids] = unique(group_prefixes, 'stable');

% Create figure with a subplot for each group
H = figure();
tiledlayout(length(unique_groups), 1, 'Padding', 'compact');

for g = 1:length(unique_groups)
    group_idx = find(group_ids == g);
    group_vars = varnames(group_idx);
    group_data = models(group_idx, :);
    group_range = prior_ranges(group_idx, :);

    % Plot
    nexttile;
    hPosterior = boxchart(group_data', 'Orientation', 'horizontal','MarkerStyle','none');
    yticklabels(group_vars);
    ax = gca;
    ax.TickLabelInterpreter = 'none';
    title(unique_groups{g});
    hold on;
    
    % Build legend dynamically
    legendHandles = hPosterior;
    legendLabels = {'posterior'};

    % Plot true values if given
    if ~isempty(p.Results.truevals)
        truevals = p.Results.truevals(group_idx);
        hTrue = plot(truevals, 1:length(group_vars), 'p', 'MarkerSize', 10, ...
             'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
        legendHandles(end+1) = hTrue;
        legendLabels{end+1} = 'true value';
    end

    % Plot best model if given
    if ~isempty(p.Results.bestmodel)
        bestmodel = p.Results.bestmodel(group_idx);
        hBest = plot(bestmodel, 1:length(group_vars), 'd', 'MarkerSize', 8, ...
             'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none');
        legendHandles(end+1) = hBest;
        legendLabels{end+1} = 'best fit';
    end

    if g == 1
        legend(legendHandles, legendLabels);
    end
    xlim([group_range(1) group_range(end)])
end

end