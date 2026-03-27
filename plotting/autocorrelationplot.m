function [H,C] = autocorrelationplot(models)
% plots autocorrelation of MCMC chains

H=figure;
[C,lags,ESS]=eacorr(models);
plot(lags,C,'.-',lags([1 end]),[0 0],'k');
grid on
xlabel('lags')
ylabel('autocorrelation');
text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
title('Markov Chain Auto Correlation')

end