function mini  = initialmodel_flatprior(prior_ranges,nWalks,varargin)
% This function determines an initial model based on the provided prior
% ranges. Initial models are drawn randomly from within the prior range.
% 
% For some models the time parameter has to decrease, so that T1 > T2 etc. 
% To make this happen enter the number of time elements that need to be 
% decreasing as scalar for VARARGIN
% Richard Ott, 2024

if nargin == 3
    dec_length  = varargin{1};
end

[nparas,~] =  size(prior_ranges);

mini = nan(nparas, nWalks);
for i = 1:nWalks
    mini(:,i) = prior_ranges(:,1) + rand(nparas,1) .* (prior_ranges(:,2) - prior_ranges(:,1));

    if nargin == 3   % make sure that time values of models decrease
        [sorted_values, ~] = sort(mini(1:dec_length,i), 'descend');
        mini(1:dec_length,i) = sorted_values;   
    end
end

end