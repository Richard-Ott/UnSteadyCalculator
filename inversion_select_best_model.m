function [best_model, best_model_like, best_walker, best_sample] = inversion_select_best_model(models, logLikeStore)
% Extract best model from logLikeStore(2, walker, sample).

posterior_like = squeeze(logLikeStore(2,:,:));
if isvector(posterior_like)
    posterior_like = reshape(posterior_like, 1, []);
end

[best_walker_like, best_sample_per_walker] = max(posterior_like, [], 2);
[best_model_like, best_walker] = max(best_walker_like);
best_sample = best_sample_per_walker(best_walker);

best_model = models(:, best_walker, best_sample);
end
