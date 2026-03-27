function N = Nforward_steadystate(E,sp,consts,Nlogical)
%NFORWARD_STEADYSTATE Calculate steady-state nuclide concentrations.
%
%   N = Nforward_steadystate(E,sp,consts)
%   N = Nforward_steadystate(E,sp,consts,Nlogical)
%
%   Computes steady-state concentrations for 10Be, 14C (and optionally 26Al)
%   for a constant erosion rate.
%
%   Inputs
%   ------
%   E : scalar or nSamp x 1 vector
%       Erosion rate in mm/ka.
%       If scalar, the same erosion rate is used for all samples.
%   sp : struct
%       Sample-specific production and attenuation parameters returned by
%       sample_parameters().
%   consts : struct
%       Constants from make_constants().
%   Nlogical : nSamp x 3 logical array
%       Mask of measured nuclides [10Be, 14C, 26Al].
%
%   Output
%   ------
%   N :
%       - column vector with entries selected by
%         Nlogical (same behavior as Nforward_discretized).
%
% Richard Ott 2026

nSamp = length(sp.P10spal);

if isscalar(E)
    E = repmat(E, nSamp, 1);
else
    E = E(:);
end


E = E ./ 1e4;  % convert mm/ka to cm/a
rho = consts.density;

att_l_10 = [consts.L_sp * ones(nSamp,1), sp.L10_nm(:), sp.L10_fm(:)];
att_l_14 = [consts.L_sp * ones(nSamp,1), sp.L14_nm(:), sp.L14_fm(:)];
att_l_26 = [consts.L_sp * ones(nSamp,1), sp.L26_nm(:), sp.L26_fm(:)];

P10 = [sp.P10spal(:), sp.P10_nm(:), sp.P10_fm(:)];
P14 = [sp.P14spal(:), sp.P14_nm(:), sp.P14_fm(:)];
P26 = [sp.P26spal(:), sp.P26_nm(:), sp.P26_fm(:)];

beta10 = rho .* E ./ att_l_10 + consts.l10;
beta14 = rho .* E ./ att_l_14 + consts.l14;
beta26 = rho .* E ./ att_l_26 + consts.l26;

N10 = sum(P10 ./ beta10, 2);
N14 = sum(P14 ./ beta14, 2);
N26 = sum(P26 ./ beta26, 2);

Nfull = [N10, N14, N26];


N = reshape(Nfull(Nlogical), [], 1);

end
