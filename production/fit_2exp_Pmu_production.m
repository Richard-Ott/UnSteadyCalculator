function [Pnmu, Lnmu, Pfmu, Lfmu] = fit_2exp_Pmu_production(Pmu,zmass,plotflag)
% Fits 2 exponentials (negative muon catpure and fast muons) to muon
% production rate profiles calculated from muon stopping rates with Balco
% functions. 
% input: Pmu - array of muon production rates with depth
%        zmass - array of depths of production rates in g/cm2
%        plotflag - do you want the fit to be plotted?
% output:
%        Pnmu - surface production rate negative muon capture (at/g/yr)
%        Lnmu - effective attenuation length negative muon capture (g/cm2)
%        Pfmu - surface production rate fast muons
%        Lfmu - effective attenuation length fast muons
%
% Richard Ott, 2025

% first-order fit to get starting point of fminsearch
pf1 = polyfit(zmass,log(Pmu),1);
P1 = exp(pf1(2));
L1 = -1./pf1(1);

Pmu = double(Pmu); % sometimes this is necessary because pressure gets passed on as single which propagates into Pmu data type

% Starting guess
options = optimset('MaxFunEvals',50000);
x0 = double([P1/1.5 P1/3 L1/2 L1.*1.5]); %[neg_P fast_P neg_L fast_L]
xopt = fminsearch(@(x) sum(((x(1).*exp(-zmass./x(3)) + x(2).*exp(-zmass./x(4)))-Pmu).^2),x0,options);

% assign output
Pnmu = xopt(1);
Pfmu = xopt(2);
Lnmu = xopt(3);
Lfmu = xopt(4);

predP = Pnmu.*exp(-zmass./Lnmu) + Pfmu.*exp(-zmass./Lfmu);
ts = 'Two exponentials';

% Do plotting

if plotflag
    figure;
    plot(Pmu,zmass,'ko','markerfacecolor','k');
    hold on;
    plot(predP,zmass,'k');
    title(ts);
    xlabel('Pmu (atoms/g/yr)');
    set(gca,'ydir','reverse','xscale','log');
    ylabel('Depth (g/cm2)');

    figure
    plot(Pmu./predP,zmass)
    set(gca,'ydir','reverse');
end

end