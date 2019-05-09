% dn_FitDN2ECoGBroadband.m


%% LOAD DATA AND PARAM FILES

dataLoc    = fullfile(dn_ECoG_RootPath, 'data');
dt_fName   = 'dn_preprocessedData.mat';
prm_fName  = 'dn_params.mat';

a = load(fullfile(dataLoc, dt_fName));
b = load(fullfile(dataLoc, prm_fName));

%% PRE-DEFINED PARAMETERS

nBoots   = 100; % bootstrap sample size
irfType  = 'uniphasic';
prm      = b.prm;

%% DERIVED PARAMETERS

bbts    = a.dt.ecog.bbts_roi;
t       = prm.ecog.t;
stim    = prm.ecog.stim;
dn      = prm.ecog.dn;

nrois   = length(bbts);

% normalize each broadband time course to its max
for iroi = 1 : nrois
   bbts{iroi} = normMax(bbts{iroi}); 
end

%% BOOTSTRAP THE ECoG DATA

bs_bbts = dn_BootstrapECoGBroadband(bbts, t, nBoots);

% save the bootstrapped data?
dn.bs_bbts = bs_bbts;

%% GRID FIT

[dn.seed, dn.seedR2] = dn_gridFit(bs_bbts', dn, stim, t, irfType);

%% FINE FIT

seed = [dn.seed, zeros(size(dn.seed, 1), 1)];
[dn.param, dn.prd, dn.r2, dn.exitflg] = dn_fineFit(bs_bbts', stim, t, dn, seed, irfType);

%% FIT THE SCALE OF THE BOOTSTRAPPED ECOG BROADBAND DATA

%% COMPUTE SUMMARY METRICS

zeroshift_param = dn.param;
zeroshift_param(:, 5) = 0; 

dn.summaryMetrics = dn_computeDerivedParams(zeroshift_param, irfType);

%% SAVE MODEL PARAMETERS

prm.ecog.dn = dn;

save(fullfile(dataLoc, prm_fName), 'prm')

%%
