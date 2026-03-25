function consts = make_constants()
%MAKE_CONSTANTS Define all constants needed by the inversion workflow.

% -------------------------------------------------------------------------
% Decay constants
% -------------------------------------------------------------------------

% Be-10 -- Chmeleff/Korschinek value, t(1/2) = 1.387 +/- 0.012 Ma
consts.l10 = -log(0.5)./1.387e6;
dldt10 = -log(0.5).*(1.387e6^-2);
consts.dell10 = sqrt((dldt10.*0.012e6)^2);

% Al-26 -- value compatible with Nishiizumi standards
% lambda = 9.83e-7 yr^-1 --> t(1/2) = 0.705 Ma; see Nishiizumi (2004)
consts.l26 = 9.83e-7;
consts.dell26 = 2.5e-8;

% C-14; t(1/2) = 5730 yr; uncertainty ~0.5 %
consts.l14 = -log(0.5)./5730;
consts.dell14 = consts.l14.*0.005;

% -------------------------------------------------------------------------
% Lal (1991) / Stone (2000) reference spallation production rates
% All values from CRONUS make_consts_v3 (Balco, BGC).
% Used by sample_parameters.m via stone2000 scaling.
% -------------------------------------------------------------------------

% Be-10 in quartz -- from recalibrating with CRONUS primary dataset 20161204
consts.P10q_St = 4.086;
consts.delP10q_St = consts.P10q_St.*0.079;  % bootstrap from secondary cal data set

% Al-26 in quartz -- from recalibrating with CRONUS primary dataset 20161204
consts.P26q_St = 28.535;
consts.delP26q_St = consts.P26q_St.*0.104;

% He-3 in quartz -- from Vermeesch data
consts.P3q_St = 115.7;
consts.delP3q_St = 10;

% He-3 in olivine/pyroxene -- from recalibrating with CRONUS primary 20160415
consts.P3op_St = 119.6;
consts.delP3op_St = consts.P3op_St.*0.11;   % bootstrap from primary data set scatter

% C-14 in quartz -- from Brent CRONUS-A data 20161204
consts.P14q_St = 12.1;
consts.delP14q_St = consts.P14q_St.*0.05;

% Ne-21 in quartz -- from SPICE data
consts.P21q_St = 16.896;
consts.delP21q_St = 1.033;

% Be-10 in pyroxene -- from Bergelin saturation data
consts.P10px_St = 3.72;
consts.delP10px_St = consts.P10px_St.*0.06; % approximately matched to Be-10 in quartz

% 8-element reference-production vector.
% Order: P3q, P3ol, P3px, P10q, P14q, P21q, P26q, P10px
% Indices 4 (P10), 5 (P14), 7 (P26) are the ones used by sample_parameters.
consts.refP_St = [ ...
	consts.P3q_St   ...   % 1 He-3 quartz
	consts.P3op_St  ...   % 2 He-3 olivine
	consts.P3op_St  ...   % 3 He-3 pyroxene
	consts.P10q_St  ...   % 4 Be-10 quartz   <-- used
	consts.P14q_St  ...   % 5 C-14  quartz   <-- used
	consts.P21q_St  ...   % 6 Ne-21 quartz
	consts.P26q_St  ...   % 7 Al-26 quartz   <-- used
	consts.P10px_St ...   % 8 Be-10 pyroxene
];

consts.delrefP_St = [ ...
	consts.delP3q_St   ...
	consts.delP3op_St  ...
	consts.delP3op_St  ...
	consts.delP10q_St  ...
	consts.delP14q_St  ...
	consts.delP21q_St  ...
	consts.delP26q_St  ...
	consts.delP10px_St ...
];

% -------------------------------------------------------------------------
% Physical constants used by Nforward_discretized
% -------------------------------------------------------------------------

consts.density = 2.65;  % g/cm3 -- typical quartz/granite density
consts.L_sp    = 160;   % g/cm2 -- spallation effective attenuation length

end