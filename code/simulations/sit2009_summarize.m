% Sit et al. 2009 result summary

% Desirable properties:
% (1) Latency:                   same across space at a given contrast, slower at lower contrast:
% (2) Rising edge:             fast in the middle (correspond to stimulus center, fovea?), slow in periphery, faster for higher contrast
% (3) Falling edge:            independent of cortical locations and stimulus contrast
% (4) Spatial distribution:  independent of contrasts

%% receptive field and stimulus parameters

sz       = 30;   % number of entries per dimension of receptive field and stimulus
rf_sig   = 1;    % the spread of the spatial receptive fielf
tau      = 0.07; % the length of the temporal receptive window
cort_sz  = 2.75; % in unit of millimeter, the extent of cortex we are looking at
stim_dur = 0.5;  % stimulus duration

% make parameters
prm  = sit2009_mkParameters(sz, rf_sig, tau, cort_sz, stim_dur);
rf   = prm.rf;
stim = prm.stim;

% create high and low contrast stimuli
contrast_levels = [1, 0.3];

hc_stim = stim; hc_stim.st = stim.st.*contrast_levels(1);
lc_stim = stim; lc_stim.st = stim.st.*contrast_levels(2);

%% Initiate computing model responses

rsp = [];

%% Compute spatio-temporal response prediction based on LINEAR and STATIC NON-LINEAR models
% This is the case that we use a static non-linear model
modelType = 'STN';

% prediction to high and low contrast stimuli
hc_rsp = sit2009_staticNonlinear(sz, hc_stim, rf, modelType);
lc_rsp = sit2009_staticNonlinear(sz, lc_stim, rf, modelType);

% linear and STN model prediction
rsp.lin = {hc_rsp.lin, lc_rsp.lin};
rsp.stn = {hc_rsp.s_nlin, lc_rsp.s_nlin};

modelType = 'ETC';
% prediction to high and low contrast stimuli
hc_rsp = sit2009_staticNonlinear(sz, hc_stim, rf, modelType);
lc_rsp = sit2009_staticNonlinear(sz, lc_stim, rf, modelType);

% ETC model prediction
rsp.etc = {hc_rsp.s_nlin, lc_rsp.s_nlin};


%% Compute spatio-temporal response prediction based on LATERAL PROPOGATION model



%% Compute spatio-temporal response prediction based on STATIC NORMALIZATION  model


%% Compute spatio-temporal response prediction based on DYNAMIC NORMALIZATION  model

modelType = 'monophasic';
hc_rsp = sit2009_DN(sz, hc_stim, rf, prm.t, modelType);
lc_rsp = sit2009_DN(sz, lc_stim, rf, prm.t, modelType);

rsp.dn_mono = {hc_rsp, lc_rsp};

modelType = 'biphasic';
hc_rsp = sit2009_DN(sz, hc_stim, rf, prm.t, modelType);
lc_rsp = sit2009_DN(sz, lc_stim, rf, prm.t, modelType);

rsp.dn_bi = {hc_rsp, lc_rsp};

%% Visualize model prediction
% 
% figure (1), clf, sit2009_visualize(rsp.lin, prm, sz, 1);
% figure (2), clf, sit2009_visualize(rsp.stn, prm, sz, 2);
% figure (3), clf, sit2009_visualize(rsp.etc, prm, sz, 3);
figure (4), clf, sit2009_visualize(rsp.dn_mono, prm, sz, 4);
figure (5), clf, sit2009_visualize(rsp.dn_bi, prm, sz, 5);