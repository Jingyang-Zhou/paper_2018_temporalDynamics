% DESCRIPTIONS

%% FIGURE SAVING KNOB

saveFigures = 0;

%% USEFUL FUNCTIONS

normMax = @(x) x./max(x);

%% PRE-DEFINED PARAMETERS

whichRois = 1 : 3; % V1 - V3;
nBoots    = 100;

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
dn_ORI  = ecogPrm.dn_ORI; % OFFSET RESPONSE INDEX
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

%% AVERAGE DATA IN EACH ECCENTRICITY BIN

ecc_bbts = {}; mecc_bbts = {};

for iroi = whichRois
    mecc_bbts{iroi} = nan(3, length(t));
    for iecc = 1 : 3
        ecc_bbts{iroi}{iecc} = normMax(bbts{iroi}(:, elec_eccIdx{iroi} == iecc));
        mecc_bbts{iroi}(iecc, :) = normMax(mean(ecc_bbts{iroi}{iecc}, 2));
    end
end

%% FIT THE BIPHASIC DN MODEL TO THE DATA

%% MAKE FOUR SETS OF SEEDS
seed = [];
% more linear, less linear, more biphasic and less biphasic

seed(1, :) = [0.02, 0.8, 0.15, 0.1, 0.05];
seed(2, :) = [0.03, 0.8, 0.1, 0.2, 0.05]; % more linear

seed(3, :) = [0.02, 0.4, 0.15, 0.1, 0.05]; % less biphasic
seed(4, :) = [0.03, 0.4, 0.1, 0.2, 0.05];

%% FIT THE BIPHASIC MODEL

param = {}; pred = {}; r2 = {};

for iroi = 1 : 3
    thisnBoots = size(bbts{iroi}, 2);
    for iset = 1 : 4
        this_seed  = repmat(seed(iset, :), [thisnBoots, 1])
        [param{iroi}(iset, :, :), pred{iroi}(iset, :, :), r2{iroi}(iset, :)] = dn_fineFit(bbts{iroi}', stim, t, bi, this_seed, irfType);
    end
    iroi
end

%% SELECT THE SET OF MODEL PREDICTION THAT PRODUCES THAT MAXIMAL R2

maxr2_idx = []; bi.param = []; bi.pred = [];

for iroi = 1 : 3
    [~, maxr2_idx{iroi}] = max(r2{iroi});
    thisnBoots = size(bbts{iroi}, 2);
    for k = 1 : thisnBoots
        bi.param{iroi}(k, :) = squeeze(param{iroi}(maxr2_idx{iroi}(k), k, :));
        bi.pred{iroi}(k, :)  = squeeze(pred{iroi}(maxr2_idx{iroi}(k), k, :));
    end
end
%% SUMMARIZE THE WEIGHT PARAMETER

mweight = []; sweight = {};

for iroi = 1 : 3
    weight = bi.param{iroi}(:, 2);
    for iecc = 1 : 3
        mweight(iroi, iecc) = mean(weight(bi.ecc.elecIdx{iroi} == iecc));
    end
end

%% PLOT PANEL A - VISUALIZE THE DIFFERENCE IN RE

fg1 = figure (1); clf
for iroi = 1 : 3
    subplot(1, 3, iroi), set(gca, 'colororder', copper(3)), hold on
    patch([0, 0.5, 0.5, 0], [0, 0, 1, 1], 0.9 * ones(1, 3));
    plot(t - 0.2, mecc_bbts{iroi}', 'linewidth', 3), 
    set(gca, 'xaxislocation', 'origin'),  axis tight, set(gca, 'fontsize', 14, 'xtick', [0, 0.5], 'ytick', [0, 1]), 
end
subplot(1, 3, 1), legend('stim.', 'low ecc.', 'mid. ecc', 'high ecc.')
fg1.Position = [1, 500, 2000, 500];

%% VISUALIZE MODEL FIT (FIGURE S4)

fgs4 =figure (100); clf
for iroi = 1 : 3
   for k = 1 : size(bbts{iroi}, 2)
      subplot(3, 15, (iroi-1) * 15 + k)
      if bi.ecc.elecIdx{iroi}(k) == 1, patch([-0.2, 1, 1, -0.2], [0, 0, 1, 1], [1, 0.9, 0.9]), % red
      elseif bi.ecc.elecIdx{iroi}(k) == 2, patch([-0.2, 1, 1, -0.2], [0, 0, 1, 1], [0.9, 0.9, 1]), % blue
      else  patch([-0.2, 1, 1, -0.2], [0, 0, 1, 1], [0.9, 1, 0.9]), % green
      end
      hold on
      patch([0, 0.5, 0.5, 0], [0, 0, 1, 1], 0.9 * ones(1, 3)), 
      plot(t - 0.2, normMax(bbts{iroi}(:, k)), 'k-', 'linewidth', 1.5),
      plot(t - 0.2, normMax(bi.pred{iroi}(k, :)), 'r-', 'linewidth', 1.5), axis tight, box off
      set(gca, 'xtick', [0, 0.5], 'ytick', [0, 1], 'xaxislocation', 'origin'), ylim([-0.2, 1])
   end
end

fgs4.Position =  [1, 500, 2700, 400];
%% PLOT PANEL B

fg2 = figure (2); clf, set(gca, 'colororder', copper(3)), hold on
patch([0, 0.5, 0.5, 0], [0, 0, 1, 1], 0.9 * ones(1, 3))

exp_prm    = [0.02, 0, 0.15, 2, 0.1, 0.05, 1];
exp_weight = [0.1, 0.75, 0.9];

for k = 1 : 3
    exp_prm(2) = exp_weight(k);
    exp_prd(k, :) = normMax(dn_DNmodel(exp_prm, stim, t));
    plot(t - 0.2, exp_prd(k, :), 'linewidth', 2.5),
end

set(gca, 'xtick', [0, 0.5], 'ytick', [0, 1], 'fontsize', 16, 'xaxislocation', 'origin'), 
legend('stim.', 'weight = 0.1', 'weight = 0.75', 'weight = 0.9'), 
xlabel('time (s)'), ylabel('model prediction'), axis tight, ylim([-0.1, 1])

fg2.Position =  [1, 250, 500, 250];

%% PLOT PANEL C

figure (3), clf, set(gca, 'colororder', gray(4)), hold on
plot(mweight', 'o', 'markersize', 10, 'linewidth', 3), set(gca, 'xtick', 1 : 3, 'xticklabel', {'low', 'mid', 'high'}, 'fontsize', 16)
xlim([0.5, 3.5]), ylim([0.4, 1]), title('weight'), axis square, xlabel('eccentricity'), ylabel('weight')

%% SAVE FIGURES

if saveFigures
    figLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
    % PANEL A
    fg1Nm = 'fg_bidnFit2ECoG_A';
    printnice(1, 0, figLoc, fg1Nm);
    % PANEL B
    fg2Nm = 'fg_bidnFit2ECoG_B';
    printnice(2, 0, figLoc, fg2Nm);
    % PANEL C
    fg3Nm = 'fg_bidnFit2ECoG_C';
    printnice(3, 0, figLoc, fg3Nm);
    % FIGURE S4
    fgs4Nm = 'fg_bidnFit2ECoG_individualFit';
    printnice(100, 0, figLoc, fgs4Nm);
end