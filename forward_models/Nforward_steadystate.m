function N = Nforward_steadystate(E,sp,consts,zm,Nlogical)
%NFORWARD_STEADYSTATE Calculate steady-state nuclide concentrations.
%
%   N = Nforward_steadystate(E,sp,consts,zm)
%   N = Nforward_steadystate(E,sp,consts,zm,Nlogical)
%
%   Computes steady-state concentrations for 10Be, 14C, and 26Al for a
%   constant erosion rate, with optional soil mixing. For zm > 0, the
%   mixed-layer solution follows the steady-state form of the box-mixing
%   equations presented by Schaller et al. (2009, Eq. 2) and Hippe (2012,
%   Appendix A, Eq. 4), adapted to the notation used in this toolbox.
%
%   Inputs
%   ------
%   E : scalar or nSamp x 1 vector
%       Erosion rate in mm/ka.
%   sp : struct
%       Sample-specific production and attenuation parameters returned by
%       sample_parameters().
%   consts : struct
%       Constants from make_constants().
%   zm : scalar
%       Soil mixing depth in cm. Use zm = 0 for the no-mixing solution.
%   Nlogical : nSamp x 3 logical array
%       Mask of measured nuclides [10Be, 14C, 26Al].
%
%   Output
%   ------
%   N : column vector with entries selected by Nlogical.
%
% Richard Ott, 2026

if nargin < 4 || isempty(zm)
    zm = 0;
end

nSamp = length(sp.P10spal);

if nargin < 5 || isempty(Nlogical)
    Nlogical = true(nSamp,3);
end

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

if zm == 0
    beta10 = rho .* E ./ att_l_10 + consts.l10;
    beta14 = rho .* E ./ att_l_14 + consts.l14;
    beta26 = rho .* E ./ att_l_26 + consts.l26;

    N10 = sum(P10 ./ beta10, 2);
    N14 = sum(P14 ./ beta14, 2);
    N26 = sum(P26 ./ beta26, 2);
else
    N10 = steady_mixed_layer(P10, att_l_10, consts.l10, rho, zm, E);
    N14 = steady_mixed_layer(P14, att_l_14, consts.l14, rho, zm, E);
    N26 = steady_mixed_layer(P26, att_l_26, consts.l26, rho, zm, E);
end

Nfull = [N10, N14, N26];
N = reshape(Nfull(Nlogical), [], 1);

end


function Nm = steady_mixed_layer(P, L, lambda, rho, Zm, E)
nSamp = size(P,1);
Nm = nan(nSamp,1);

for s = 1:nSamp
    beta = lambda + (rho .* E(s)) ./ L(s,:);
    exp_fac = exp(-(rho .* Zm) ./ L(s,:));

    Pbar = sum(P(s,:) .* (L(s,:) ./ (rho .* Zm)) .* (1 - exp_fac));
    Nb_zm = sum(exp_fac .* (P(s,:) ./ beta));

    r = lambda + E(s) ./ Zm;
    Nm(s) = (Pbar + (E(s) ./ Zm) .* Nb_zm) ./ r;
end
end
