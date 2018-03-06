% dn_fitDNbi2ECoG_v2

%% PRE-DEFINED PARAMETERS

whichRois = 1 : 3; % V1 - V3;
figureLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
nBoots    = 100;

%% USEFUL FUNCTIONS

normMax = @(x) x./max(x);

%% LOAD DATA AND PARAM FILES

dataLoc    = fullfile(dn_ECoG_RootPath, 'data');
dt_fName   = 'dn_preprocessedData.mat';
prm_fName  = 'dn_params.mat';

a = load(fullfile(dataLoc, dt_fName));
b = load(fullfile(dataLoc, prm_fName));

%% DERIVED PARAMETERS

ecog    = a.dt.ecog;
ecogPrm = b.prm.ecog;

ElecIdx = ecog.elec_roi;
ecc     = ecog.prf.ecc;
bbts    = ecog.bbts_roi;
nrois   = length(ecc);

t       = ecogPrm.t;
stim    = ecogPrm.stim;
dn      = ecogPrm.dn;
dn_ORI  = ecogPrm.dn_ORI;
bi      = ecogPrm.dn_bi;

%% MAKE ECCENTRICITY BINS

bin(1, :) = [0, 5]; bin(2, :) = [5, 10]; bin(3, :) = [10, inf];

%% ASSIGN ELECTRODES TO ECCENTRICITY BINS

elec_eccIdx = {};

% COMPUTE WHICH ELECTROD CORRESPONDS TO WHICH ECCENTRICITY BIN ------------
for iroi = whichRois
    idxLow = ecc{iroi} < 5;    elec_eccIdx{iroi}(idxLow) = 1;
    idxMid = (ecc{iroi} >= 5) & (ecc{iroi} < 10); elec_eccIdx{iroi}(idxMid) = 2;
    idxHigh = ecc{iroi} >= 10; elec_eccIdx{iroi}(idxHigh) = 3;
end
% SAVE THIS PART OF THE ANALYSIS ------------------------------------------
bi.ecc.bin     = bin;
bi.ecc.elecIdx = elec_eccIdx;

%% GROUP ELECTRODES AND THEIR TIME COURSES ACCORDINGLY, AND BOOTSTRAP

ecc_bbts = {}; ecc_bbtsBs = {};

for iroi = whichRois
    for iecc = 1 : 3
        ecc_bbts{iroi}{iecc} = normMax(bbts{iroi}(:, elec_eccIdx{iroi} == iecc));
    end
     bs_bbtsBs{iroi} = dn_BootstrapECoGBroadband(ecc_bbts{iroi}, t, nBoots);
end

%% VISUALIZE THE AVERAGE TIME COURSE (AVERAGED ACROSS ELECTRODES WITHIN EACH ECCENTRICITY BIN)

figure (1), clf
for iroi = whichRois
    subplot(1, 3, iroi), set(gca, 'colororder', copper(3)), hold on
    % PLOT THE STIMULUS -----------------------------------------------
    patch([0, 0.5, 0.5, 0], [0, 0, 1, 1], 0.9 * ones(1, 3));
    for iecc = 3 : -1 : 1
        idx = (iecc - 1) * nBoots + 1 : iecc * nBoots;
        m = normMax(mean(bs_bbtsBs{iroi}(:, idx), 2));
        s = std(bs_bbtsBs{iroi}(:, idx), [], 2);
        plot(t-0.2, m, 'linewidth', 3),
    end
    set(gca, 'xaxislocation', 'origin'),  axis tight, set(gca, 'fontsize', 14, 'xtick', [0, 0.5], 'ytick', [0, 1])
end
subplot(1, 3, 1), legend('stim.', 'high ecc.', 'mid. ecc', 'low ecc.')

%% FIT THE BIPHASIC MODEL TO THE BOOTSTRAPPED TIME COURSES

% RESHAPE THE TIME COURSES
bs_bbtsBs{3} = bs_bbtsBs{3}(:, [1 : 100, 201 : 300]);

% GRID FIT THE BIPHASIC DN MODEL 
irfType = 'biphasic';
for iroi = whichRois
    [bi.seed{iroi}, bi.seedR2{iroi}] = dn_gridFit(bs_bbtsBs{iroi}', bi, stim, t, irfType);
end

%% FINE FIT THE BIPHASIC MODEL TO THE ECC. DATA

bi.param = {}; bi.pred = {}; bi.r2 = {}; bi.exitFlg = {};

for iroi = whichRois
    thisnBoots  = size(bs_bbtsBs{iroi}, 2);
    seed = repmat([0.05, 0, 0.1, 0.1, 0], [thisnBoots, 1]);
    seed(:, 2) = bi.seed{iroi};
    [bi.param{iroi}, bi.pred{iroi}, bi.r2{iroi}, bi.exitFlg{iroi}] = dn_fineFit(bs_bbtsBs{iroi}', stim, t, bi, seed, irfType);
    iroi
end

%% VISUALIZE THE MODEL FIT
