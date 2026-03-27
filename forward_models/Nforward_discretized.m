function N = Nforward_discretized(E,T,sp,consts,scenario,Nlogical,varargin)
% This function calculates concentrations N10 and N14 for mutiple step 
% changes in erosion and for multiple samples at once.
%
% Input:
%       - E: erosion rates in mm/ka from old to current for n samples for 
%           m step changes. Input format:
%               - n x m for type 'step'
%               - 1 x n for type 'samestep'
%               - scalar for type 'spike' (background erosion rate)
%        
%       - T: 1 x m vector of step change timings (yrs BP)
%           from old to young
%       - sp: sample specific parameters (Pspal & pressure)
%       - consts: the consts_v3.mat file from Cronus v3
%       - precalculated muon production rates from Cronus v3 (Nmu.mat)
%       - type: type of erosion scenario. Either: 
%               - 'step': step changes in erosion rates with individual
%               rates per sample
%               - 'samestep': step changes in erosion with common change
%               factor
%               - 'spike': spikes of soil loss
%       - varargin: empty for 'step'
%
% Output:
%       - [N10,N14]: concentrations at/g of 10Be and 14C
%
% Richard Ott & Dirk Scherler, 2024

nSamp = length(sp.P10spal);

E = E./1e4;   % convert to cm/a
T_time_spans   = [inf, diff(T')*(-1)];        % time span of every time interval between erosion changes

%% reshape erosion rate array and calculate depths at steps for different scenarios

switch scenario
    case 'step'
        segment_depths = E.*T_time_spans;            % the exhumation occuring during every erosion time interval in cm

    case 'samestep'
        CHG = varargin{1};
        % make matrix of erosion rates for different time steps by multiplying E
        % with the change factors
        E = [E, repmat(E,1,length(CHG))];
        E(:,2:end) = E(:,2:end) * diag(CHG);
        segment_depths = E.*T_time_spans;            % the exhumation occuring during every erosion time interval in cm

    case 'samebackground_step'
        CHG = varargin{1};
        E = repmat(E,nSamp,length(T));
        E(:,2:end) = E(:,2:end) .* CHG;
        segment_depths = E.*T_time_spans;            % the exhumation occuring during every erosion time interval in cm

    case 'samebackground_samestep'
        CHG = varargin{1};
        E = repmat(E,nSamp,length(T));
        E(:,2:end) = E(:,2:end) * diag(CHG);
        segment_depths = E.*T_time_spans;            % the exhumation occuring during every erosion time interval in cm

    case 'spike'
        loss = varargin{1};
        E = repmat(E,1, length(T));                      % make matrix of erosion rates for easy calling in concentration loop
        segment_depths = E.*T_time_spans+[inf(nSamp,1) loss];     % calcuate the exhumation occuring during every erosion time interval in cm

    case 'samespike'
        loss = varargin{1};
        E = repmat(E,1, length(T));
        segment_depths = E.*T_time_spans+[inf(nSamp,1) repmat(loss',nSamp,1)];     % calcuate the exhumation occuring during every erosion time interval in cm

    case 'samebackground_spike'
        loss = varargin{1};
        E = repmat(E,nSamp,length(T));                      % make matrix of erosion rates for easy calling in concentration loop
        segment_depths = E.*T_time_spans+[inf(nSamp,1) loss];

    case 'samebackground_samespike'
        loss = varargin{1};
        E = repmat(E,nSamp,length(T));                      % make matrix of erosion rates for easy calling in concentration loop
        segment_depths = E.*T_time_spans+[inf(nSamp,1) repmat(loss',nSamp,1)];
    case 'curve'
        scaleFactor = varargin{1};

        absolute_erosion_change = sp.curvechange .* scaleFactor; % calculate erosion rate change from relative differences from curve and scaling factor
        % make matrix of erosion rates for different time steps by multiplying E with the change factors
        E = repmat(E,1,length(T));
        E(:,2:end) = E(:,2:end)+ E(:,2:end) .* absolute_erosion_change;
        segment_depths = E.*T_time_spans;            % the exhumation occuring during every erosion time interval in cm
end

t_depths       = [fliplr(cumsum(fliplr(segment_depths(:,2:end)),2)), zeros(nSamp,1)];     % depths at every step change T in cm


%% start cosmo calulcation

% assign variables --------------------------------------------------------
rho = consts.density;  

% attentuation lengths
att_l_10(:,1) = consts.L_sp* ones(nSamp,1);          % spallation
att_l_10(:,2) = sp.L10_nm;        % negative muon 
att_l_10(:,3) = sp.L10_fm;        % fast muon 
att_l_14(:,1) = consts.L_sp* ones(nSamp,1);          % spallation
att_l_14(:,2) = sp.L14_nm;        % negative muon 
att_l_14(:,3) = sp.L14_fm;        % fast muon 
att_l_26(:,1) = consts.L_sp* ones(nSamp,1);          % spallation
att_l_26(:,2) = sp.L26_nm;        % negative muon 
att_l_26(:,3) = sp.L26_fm;        % fast muon 

% production rates 
P10(:,1) = sp.P10spal; 
P14(:,1) = sp.P14spal; 
P26(:,1) = sp.P26spal; 
P10(:,2) = sp.P10_nm; 
P14(:,2) = sp.P14_nm; 
P26(:,2) = sp.P26_nm; 
P10(:,3) = sp.P10_fm; 
P14(:,3) = sp.P14_fm; 
P26(:,3) = sp.P26_fm; 

% do we need to calculate Alumimium (currently I run 26Al as a separate,
% assuming that it will be the least measured nuclide, so instead of
% modyfying the current loop, I just add a new loop with an if statement)
% taking Al into the main loop would be cleaner, but I'm trying to optimize
% for speed in the inversion.
Al = any(Nlogical(:,3));

% calculate concentrations --------------------------------------------
N10 = 0;
N14 = 0; 

% loop through production pathways (first spallation, then negative and fast muons)
for i = 1:3 
    
    % segments of production profile
    N10i = 0;
    N14i = 0;
    
    for j = 1:length(T)  % loop through erosion segments

        beta10 = rho .* E(:,j) ./ att_l_10(:,i) + consts.l10;
        beta14 = rho .* E(:,j) ./ att_l_14(:,i) + consts.l14;
    
        N10_segment = P10(:,i)./beta10 .* ...                    % production
            exp(-(rho .* t_depths(:,j) ./ att_l_10(:,i))) .* ...    % depth cut-off
            (1 - exp(-beta10.*T_time_spans(j)));              % time to develop the concentration profile
        N14_segment = P14(:,i)./beta14 .* ...
            exp(-(rho .* t_depths(:,j) ./ att_l_14(:,i))) .* ...
            (1 - exp(-beta14.*T_time_spans(j)));

        N10i = N10i .* exp(-consts.l10.*T_time_spans(j)) + N10_segment; % add segments and don't forget decay
        N14i = N14i .* exp(-consts.l14.*T_time_spans(j)) + N14_segment; % add segments and don't forget decay       
    end

    N10 = N10 + N10i;
    N14 = N14 + N14i;
end
 

if Al       % calculate Al separetely assuming that it will be used the least and we run code faster if we dont calculate this
N26=0;
    for i = 1:3 
    N26i = 0;    
    for j = 1:length(T)  % loop through erosion segments
        beta26 = rho .* E(:,j) ./ att_l_26(:,i) + consts.l26;

        N26_segment = P26(:,i)./beta26 .* ...                    % production
            exp(-(rho .* t_depths(:,j) ./ att_l_26(:,i))) .* ...    % depth cut-off
            (1 - exp(-beta26.*T_time_spans(j)));              % time to develop the concentration profile

        N26i = N26i .* exp(-consts.l26.*T_time_spans(j)) + N26_segment; % add segments and don't forget decay
    end
    N26 = N26 + N26i;
    end
end

if Al
    Ncalculated = [N10;N14;N26];
else
    Ncalculated = [N10;N14; zeros(nSamp,1)];
end
N = Ncalculated(Nlogical);  % only pass on nuclide values that were measured
end
