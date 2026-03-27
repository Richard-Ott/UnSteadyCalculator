function logical_prior = flatprior_multistep(model, prior_ranges)
% this function outputs the log-prior for the MCMC given a
% multi-step chnage in erosion rates at times T. 
% This function is for flat prior. Therefore it outputs a logical (0 - out
% of prior range, 1 - within prior range).
% Richard Ott, 2024


logical_prior = and(model > prior_ranges(:,1), model < prior_ranges(:,2));

end