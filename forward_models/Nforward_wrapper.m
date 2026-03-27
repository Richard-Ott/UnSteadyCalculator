function N = Nforward_wrapper(model,sp,consts,scenario,nsteps,Nlogical)
% This is a wrapper function for Nforward_discretized for easy use with the
% MCMC algorithm.
% Richard Ott, 2024
%
% Input:
%       - model: inversion model as column vector
%       - sp: sample specific parameters (Pspal & pressure)
%       - consts: the consts_v3.mat file from Cronus v3
%       - precalculated muon production rates from Cronus v3 (Nmu.mat)
%       - scenario: type of erosion scenario. 
%       - nsteps:  the number of step changes or spike erosion losses, = 0
%       for curve model
%       - Nlogical: a logical table n x 3 that shows which nuclides were
%       measured for which sample (column order: 10Be, 14C, 26Al)
%
% Output:
%       - N = [N10,N14]: concentrations at/g of 10Be and 14C
%
% Richard Ott 2024

nSamp = length(sp.P10spal);

T = [model(1:nsteps); 0];
if ~strcmp(scenario,'step')
    E = model(nsteps+1:end-nsteps);
    changevar = model(end-nsteps+1:(end));
end
%% reshape input
switch scenario
    case 'step'
        E = model(nsteps+1:end);
        E = reshape(E,[nSamp, length(T)]);

    case 'samestep'
        E = model(nsteps+1:end-nsteps);
        changevar = model(end-nsteps+1:(end));

    case 'samebackground_step'
        E = model(nsteps+1);
        changevar = model(nsteps+2:end);
        changevar = reshape(changevar,[nSamp, length(T)-1]);

    case 'samebackground_samestep'
        E = model(nsteps+1);
        changevar = model(nsteps+2:end);

    case 'spike'
        E = model(nsteps+1: nsteps+nSamp);
        changevar = model(nsteps+nSamp+1:end);
        changevar = reshape(changevar,[nSamp, length(T)-1]);

    case 'samespike'
        E = model(nsteps+1: nsteps+nSamp);
        changevar = model(nsteps+nSamp+1:end);

    case 'samebackground_spike'
        E = model(nsteps+1);
        changevar = model(nsteps+2:end);
        changevar = reshape(changevar,[nSamp, length(T)-1]);

    case 'samebackground_samespike'
        E = model(nsteps+1);
        changevar = model(nsteps+2:end);
    case 'curve'
        if length(model) > 2*nSamp    % initial models can include time values while sampler models do not
            E = model(nsteps+1:nsteps+nSamp);
            changevar = model(nsteps+nSamp+1:end);
        else
            E = model(1:+nSamp);
            changevar = model(nSamp+1:end);
        end
        T = [sp.t'; 0];
end

%% run forward model

if strcmp(scenario,'step')
    N = Nforward_discretized(E,T,sp,consts,scenario,Nlogical);
else
    N = Nforward_discretized(E,T,sp,consts,scenario,Nlogical,changevar);
end


end
