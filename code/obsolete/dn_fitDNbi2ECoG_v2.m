% dn_fitDNbi2ECoG_v2
%
% DESCRIPTIONS ------------------------------------------------------------
% In this file, I am going to make all the analysis for figure 4 (the
% eccentricity analysis figure).
% Here are the analyses steps: 
% (1) Find eccentricity bins, and identify which bin each electrode belongs
%     to.
% (2) For the first panel, I am going to show eccentricity bin-averaged
%     time courses for V1 - V3, and the correpsonding computed offset
%     transient index.
% (3) For the second panel, I am going to fit the biphasic DN model to each
%     electrodes, and show the correpsonding weight parameter for each bin.


%% PRE-DEFINED PARAMETERS

whichRois = 1 : 3; % V1 - V3;
figureLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
nBoots    = 100;

%% USEFUL FUNCTIONS


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
ecc     = ecog.ecc;
bbts    = ecog.bbts_roi;
nrois   = length(ecc);

t       = ecogPrm.t;
stim    = ecogPrm.stim;
dn      = ecogPrm.dn;
dn_OTI  = ecogPrm.dn_ORI;
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

%% GROUP ELECTRODES AND THEIR TIME COURSES ACCORDINGLY

bi.ecc_bbts = [];

for iroi = whichRois
    this_roi = normMax(bbts{iroi});
    for iecc = 1 : 3
        bi.ecc_bbts(iroi,iecc, :) = mean(this_roi(:, elec_eccIdx{iroi} == iecc), 2);
    end
end

%% VISUALIZE ECCENTRICITY BINNED ROI TIME SERIES

mecc_bbts = normMax(squeeze(nanmean(bi.ecc_bbts))');

figure (1), clf, set(gca, 'colororder', copper(3)), hold on
patch([0, 0.5, 0.5, 0], [0, 0, 1, 1], 0.8 * ones(1, 3))
patch([0.5, 0.8, 0.8, 0.5], [0, 0, 1, 1], [0.8, 1, 0.8])
plot(t-0.2, mecc_bbts', 'linewidth', 3), axis tight,
set(gca, 'xaxislocation', 'origin'), set(gca, 'fontsize', 16 ,'xtick', [0, 0.5, 0.8], 'ytick', [0, 1])

figLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
fig1Nm = 'ana_eccBinBroadbandAveraged';

printnice(1, 0, figLoc,fig1Nm);

% MAKE A FIGURE SEPARATE FOR EACH ROI -------------------------------------
title_txt = {'V1', 'V2', 'V3'};

figure (2), clf
for iroi = whichRois
   toplot = normMax(squeeze(bi.ecc_bbts(iroi, :, :))');
   subplot(1, 3, iroi), set(gca, 'colororder', copper(3)), hold on
   patch([0, 0.5, 0.5, 0], [0, 0, 1, 1], 0.8 * ones(1, 3))
   patch([0.5, 0.8, 0.8, 0.5], [0, 0, 1, 1], [0.8, 1, 0.8])
   plot(t-0.2, toplot, '-', 'linewidth', 3), axis tight, ylim = ([-0.1, 1]), axis square
   set(gca, 'xaxislocation', 'origin'), set(gca, 'fontsize', 16 ,'xtick', [0, 0.5], 'ytick', [0, 1])
   title(title_txt{iroi})
end
figLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
fig2Nm = 'ana_eccBinBroadbandPerRoi';

printnice(2, 0, figLoc,fig2Nm);

%% COMPUTE OFFSET TRANSIENT INDEX
% 
% OTI = [];
% 
% x0 = 0.01;
% 
% % Fit an exponential decay to [0.2 0.7] and [1, 1.2] time range of the
% % data. Then compute the difference between the data and the model
% % prediction.
% fit_rng = [201 : 700, 1000 : 1204];
% 
% for iroi = whichRois
%     for iecc = 1 : 3
%         tofit = normMax(squeeze(bi.ecc_bbts(iroi, iecc, :))');
%         if ~isnan(tofit),
%             % find maximum
%             maxIdx(iroi, iecc) = find(tofit == max(tofit), 1);
%             fit_rng = [maxIdx(iroi, iecc) : 700, 1000 : 1204];
%             OTI.prm(iroi, iecc) = fminsearch(@(x) dn_fitExponential(x, tofit, fit_rng, t), x0);
%             OTI.prd{iroi, iecc} = exp(-t(maxIdx(iroi, iecc) : end)./OTI.prm(iroi, iecc));
%             OTI.prd{iroi, iecc} = normMax(OTI.prd{iroi, iecc});
%         end
%     end
% end

%% VISUALIZE THE EXPONENTIAL PREDICTION

% figure (3), clf
% 
% for iroi = whichRois
%     subplot(1, 3, iroi)
%     for iecc = 1 : 3
%         dt2plot = normMax(squeeze(bi.ecc_bbts(iroi, iecc, :)));
%         plot(t, dt2plot + iecc, 'k-', 'linewidth', 2), hold on
%         % PLOT EXPONENTIAL  PREDICTION
%         if ~isempty(OTI.prd{iroi, iecc})
%             thismax = maxIdx(iroi, iecc);
%             rng = thismax : length(t);
%             plot(t(rng), iecc + OTI.prd{iroi, iecc}, 'r', 'linewidth', 2)
%         end
%     end
%     axis tight, box off
% end

%% THE SECOND WAY TO COMPUTE THE OFFSET RESPONSE INDEX: FIT THE DN MODEL TO THE WHOLE TIME COURSE
trimIdx = 1 : 700; stimTrim = stim(trimIdx); tTrim = t(trimIdx); irfType = 'uniphasic';

for iroi = whichRois
    nElec = length(ecc{iroi});
    toFit = bbts{iroi}';
    [ORI.seed{iroi}, ORI.seedR2{iroi}] = dn_gridFit(toFit, dn, stim, t, irfType);
end

% FINE FIT ----------------------------------------------------------------
for iroi = whichRois
    nElec = length(ecc{iroi});
    toFit = bbts{iroi}';
    seed  = [ORI.seed{iroi}, zeros(nElec, 1)];
    seedTrim = [ORI.seedTrim{iroi}, zeros(nElec, 1)];
    
    % COMPUTE MODEL FIT TO BOTH THE FIRST 700MS AND THE ENTIRE TIME COURSE
   [param_fl{iroi}, ~, ~, exitflg_fl{iroi}] = dn_fineFit(toFit, stim, t, dn_OTI, seed, irfType);
 
end

%% FIT THE BIPHASIC MODEL TO THE BOOTSTRAPPED TIME COURSES

% RESHAPE THE TIME COURSES
bs_bbtsBs{3} = bs_bbtsBs{3}(:, [1 : 100, 201 : 300]);

irfType = 'biphasic';

for iroi = whichRois
    [bi.seed{iroi}, bi.seedR2{iroi}] = dn_gridFit(bs_bbtsBs{iroi}', bi, stim, t, irfType);
    
end

%% FINE FIT
bi.param = {}; bi.pred = {}; bi.r2 = {}; bi.exitFlg = {};

for iroi = whichRois
    thisnBoots  = size(bs_bbtsBs{iroi}, 2);
    seed = repmat([0.05, 0, 0.1, 0.1, 0], [thisnBoots, 1]);
    seed(:, 2) = bi.seed{iroi};
    [bi.param{iroi}, bi.pred{iroi}, bi.r2{iroi}, bi.exitFlg{iroi}] = dn_fineFit(bs_bbtsBs{iroi}', stim, t, bi, seed, irfType);
end

%% VISUALIZE THE MODEL FIT
