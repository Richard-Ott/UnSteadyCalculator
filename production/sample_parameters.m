function [sp] = sample_parameters(lat,lon,altitude,consts)
% get top spallation production from Cronus v3 for 10Be, 14C, 26Al
% multi sample input possible (provide variables as vectors)
% calculates muon surface production rates and effective attenuation
% lengths using Blaco 2017 model 1A code.
%
% Richard Ott, 2024

% addpath('C:\Users\rott\OneDrive - UvA\Richard\PhD_ETH\matlab\CRONUS cosmo calculation\cronus 3.0\online-calculators-v3')
% addpath('C:\Users\rott\OneDrive - UvA\Richard\PhD_ETH\matlab\CRONUS cosmo calculation\cronus 3.0\online-calculators-v3\data')

sp.pressure = ERA40atm(lat, lon, altitude);  % atmospheric pressure

if isrow(sp.pressure)   % make sure its a column vector
    sp.pressure = sp.pressure';
end

P10_SLHL = consts.refP_St(4);                % Stone scaling ref spallation production
P14_SLHL = consts.refP_St(5);
P26_SLHL = consts.refP_St(7);

sf_St = (stone2000(lat,sp.pressure,1))';        % scaling factor

sp.P10spal = sf_St .* P10_SLHL;              % site specific spallation production
sp.P14spal = sf_St .* P14_SLHL;
sp.P26spal = sf_St .* P26_SLHL;


%% Balco 2017 muon production rates

% Be-10, model 1A
const10A.Natoms = 2.006e22;   % oxygen atoms per g quartz
const10A.k_neg = 0.00191 .* 0.704 .* 0.1828; % From BCO fit
const10A.sigma0 = 0.280e-30; % From BCO fit

% C-14, model 1A
fstar14 = 0.116; % Balco 2017 Model 1A fit. Lupker 2015 finds 0.124, Heisinger 2022 0.139
const14A.Natoms = 2.006e22;
const14A.k_neg = fstar14 .* 0.704 .* 0.1828; % From BCO fit,   f_star*f_C*f_D
const14A.sigma0 = 0.45e-27./190; % Use Heisinger value, for alpha = 1

% Al-26, model 1A
const26A.Natoms = 1.003e22;
const26A.k_neg = 0.0133 .* 0.296 .* 0.6559; % From BCO fit
const26A.sigma0 = 3.89e-30; % From BCO fit

sp.D = 25;   % depth of profile (m)
sp.Nn = 100; % number of nodes
z = linspace(0,1,sp.Nn);
sp.z = sp.D*z.^2; % depth node spacing
sp.zmass = sp.z*100 * consts.density;  % convert depth to g/cm2 for Balco code

% calculate muon production depth profiles
for i = 1:length(sp.pressure)
    sp.P10_mu = P_mu_total_alpha1(sp.zmass,sp.pressure(i),const10A,'no');
    sp.P14_mu = P_mu_total_alpha1(sp.zmass,sp.pressure(i),const14A,'no');
    sp.P26_mu = P_mu_total_alpha1(sp.zmass,sp.pressure(i),const26A,'no');

    % fit expontentials
    [sp.P10_nm(i), sp.L10_nm(i), sp.P10_fm(i), sp.L10_fm(i)] = ...
        fit_2exp_Pmu_production(sp.P10_mu,sp.zmass,false);
    [sp.P14_nm(i), sp.L14_nm(i), sp.P14_fm(i), sp.L14_fm(i)] = ...
        fit_2exp_Pmu_production(sp.P14_mu,sp.zmass,false);
    [sp.P26_nm(i), sp.L26_nm(i), sp.P26_fm(i), sp.L26_fm(i)] = ...
        fit_2exp_Pmu_production(sp.P26_mu,sp.zmass,false);

end

end