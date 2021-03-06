%% dn_mkParam

%% SAVE KNOB

saveData = 1;

%% LOAD PARAMS FILE

dataLoc    = fullfile(dn_ECoG_RootPath, 'data');
prm_fName  = 'dn_params.mat';

b = load(fullfile(dataLoc, prm_fName));

%% PRE-DEFINED VARIABLES

nSteps = 10;

%% DERIVED PARAMETERS

prm = b.prm;

%% MAKE ECOG DN (MODEL FIT) PARAMETERS

ecog = prm.ecog;

% FOR GRID FIT ------------------------------------------------------------
ecog.dn.tau1steps  = linspace(0.07, 0.5, nSteps);
ecog.dn.tau2steps  = linspace(0.07, 0.5, nSteps);
ecog.dn.nsteps     = linspace(1, 6, nSteps);
ecog.dn.sigmasteps = linspace(0.01, 0.2, nSteps);

% FOR FINE FIT ------------------------------------------------------------
ecog.dn.fitOptions = optimset('Display', 'notify', 'Algorithm', 'trust-region-reflective', 'MaxIter', 2000);
ecog.dn.lb         = [0.07, 0.07, 1, 0.01, 0.001];
ecog.dn.ub         = [1, 1, 6, 0.5, 0.1];

% RE-SAVE TO THE ORIGINAL STRUCTURE ---------------------------------------
prm.ecog = ecog;

%% MAKE ECOG DNbi (MODEL FIT PARAMETERS)

% FOR GRID FIT ------------------------------------------------------------
ecog = prm.ecog;
ecog.dn_bi.weightsteps = linspace(0, 1, nSteps);
ecog.dn_bi.lb = [0, 0, 0, 0.001, 0]; % tau1, weight, tau2, sigma, shift
ecog.dn_bi.ub = [0.5, 1, 0.5, 1, 0.5];

% RE-SAVE TO THE ORIGINAL STRUCTURE ---------------------------------------
prm.ecog = ecog;

%% FOR COMPUTING OFFSET RESPONSE INDEX ------------------------------------

ecog.dn_ORI.lb = [0.05, 0.05, 0, 0.001, 0];
ecog.dn_ORI.ub = [0.3, 0.5, 6, 1, 0.01];

% RE-SAVE TO THE ORIGINAL STRUCTURE ---------------------------------------
prm.ecog = ecog;

%% ECoG MEASUREMENT PARAMETERS

%(The following parameters are made in file dn_chooseElectrodes.m)

% prm.ecog:

% NAMES OF ALL THE ROIs (COMBINED ACROSS VISUAL FIELD MAPS) ---------------
% prm.ecog.roiNm = {'V1', 'V2', 'V3', 'lateral', 'ventral', 'dorsal'};
%                In the text, we also use "LA", "VA" and "DA" to represent lateral,
%                ventral, and dorsal areas.
%
% ELECTRODE SELECTION CRITERION -------------------------------------------
% prm.ecog.selThresh = 1;
%                 We choose the electrodes the maximum response of which
%                 is over at least 100% increase compared to the baseline
%
% prm.ecog.tBase     = 1 : 200;
%                 The range of time points that is considered baseline.

%% 
if isfield(prm, 'mua'), mua = prm.mua;
else prm.mua = [];
end
mua.dn.tau1steps  = linspace(0.02, 0.3, nSteps);
mua.dn.tau2steps  = linspace(0.02, 0.3, nSteps);
mua.dn.nsteps     = linspace(1, 6, nSteps);
mua.dn.sigmasteps = linspace(0.01, 0.2, nSteps);

mua.dn.lb         = [0.02, 0.02, 1, 0.001, 0.001];
mua.dn.ub         = [1, 1, 6, 1, 0.01];

prm.mua = mua;

%% SAVE DATA

if saveData,
   save(fullfile(dataLoc, prm_fName), 'prm') 
end
