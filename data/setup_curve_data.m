function [sp, cfg] = setup_curve_data(sp,cfg)
% loads pollen curve data for Crete for the inversion
% apply pollen data as test

curvedata = load('./data/pollen.mat');
pollen = curvedata.pollen;

timebreaks = [10000, 6200, 700, 0];
for i = 1:length(timebreaks) - 1
    timeRange = pollen.yearsBP >= timebreaks(i+1) & pollen.yearsBP < timebreaks(i); % Define the time range for this period
    meanPercTree(i) = mean(pollen.percTree(timeRange)); % Calculate the mean percTree for this time range
end
noTreePerc = 100-meanPercTree; 
curvechanges    = (noTreePerc(2:end) ./ noTreePerc(1) -1);

sp.t = timebreaks(2:end-1);
sp.curvechange = curvechanges;            % these are the base changes that will be scaled later on

cfg.nsteps = length(sp.t);   % adjust number of time steps

end