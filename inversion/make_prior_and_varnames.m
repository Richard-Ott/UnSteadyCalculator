function [prior_range,var_names] = make_prior_and_varnames(scenario,T,E,LOSS,CHG,n,steps)
% This function makes the maktrix of prior ranges and the variable names
% used later. Currently, all priors of a certain type (time, erosion, loss,
% changefactor) have the same prior range. But one could adapt this
% function to use different prior ranges for different, e.g., erosion
% steps.
% Richard Ott, 2024

time_prior = [repmat(T,steps,1)];
var_names  = arrayfun(@(x) sprintf('T%d', x), 1:steps, 'UniformOutput', false); % time var names

switch scenario
    case 'step'
        prior_range = [time_prior; repmat(E, n*(steps+1), 1)];
        for i = 1:steps+1
            var_names = [var_names, arrayfun(@(x) sprintf('E%d_sample%d', i, x), 1:n, 'UniformOutput', false)];
        end
    case 'samestep'
        prior_range = [time_prior; repmat(E, n, 1); repmat(CHG,steps,1)];
        var_names = [var_names, arrayfun(@(x) sprintf('E1sample%d', x), 1:n, 'UniformOutput', false),...
        arrayfun(@(x) sprintf('ChangeFactor%d', x), 1:steps, 'UniformOutput', false)];
    case 'samebackground_step'
        prior_range = [time_prior; E; repmat(CHG,n,1)];  
        var_names = [var_names, 'E', arrayfun(@(x) sprintf('ChangeFactor%d', x), 1:n, 'UniformOutput', false)];
    case 'samebackground_samestep'
        prior_range = [time_prior; E; CHG];  
        var_names = [var_names, 'E','ChangeFactor'];
    case 'spike'
        prior_range = [time_prior; repmat(E, n, 1); repmat(LOSS, n*steps, 1)];
        var_names = [var_names, arrayfun(@(x) sprintf('E_sample%d', x), 1:n, 'UniformOutput', false)];
        for i = 1:steps
            var_names = [var_names, arrayfun(@(x) sprintf('Loss%d_sample%d', i, x), 1:n, 'UniformOutput', false)];
        end
    case 'samespike'
        prior_range = [time_prior; repmat(E, n, 1); repmat(LOSS,steps,1)];
        var_names = [var_names, arrayfun(@(x) sprintf('Esample%d', x), 1:n, 'UniformOutput', false),...
        arrayfun(@(x) sprintf('Loss%d', x), 1:steps, 'UniformOutput', false)];
    case 'samebackground_spike'
        prior_range = [time_prior; E; repmat(LOSS,n,1)];  
        var_names = [var_names, 'E', arrayfun(@(x) sprintf('Loss%d', x), 1:n, 'UniformOutput', false)];
    case 'samebackground_samespike'
        prior_range = [time_prior; E; repmat(LOSS,steps,1)];  
        var_names = [var_names, 'E', arrayfun(@(x) sprintf('Loss%d', x), 1:steps, 'UniformOutput', false)];
    case 'curve'
        prior_range = [repmat(E, n, 1); repmat(CHG,n,1)];
        var_names = [var_names,arrayfun(@(x) sprintf('E1sample%d', x), 1:n, 'UniformOutput', false),...
        arrayfun(@(x) sprintf('ChangeFactor%d', x), 1:n, 'UniformOutput', false)];
end

end