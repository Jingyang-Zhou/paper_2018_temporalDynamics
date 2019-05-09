
% dn_FitDN_bi3ECoGEcc

%% PRE-DEFINED PARAMETERS

whichRois = 1 : 3; % V1 - V3;
figureLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');

%% USEFUL FUNCTIONS

normMax = @(x) x./max(x);

%% LOAD DATA AND PARAM FILES

dataLoc    = fullfile(dn_ECoG_RootPath, 'data');
dt_fName   = 'dn_preprocessedData.mat';
prm_fName  = 'dn_params.mat';

a = load(fullfile(dataLoc, dt_fName));
b = load(fullfile(dataLoc, prm_fName));

%% MAKE ELECTRODE ECC

% c = load('/Volumes/server/Projects/ECoG/Figures/pRF_data/Power/oneBar/CAR/Subj17_exp_bb.mat');
% x0 = c.params.params(:, 1);
% y0 = c.params.params(:, 2);
% s0 = c.params.params(:, 3);
%
% subj = 17; nrois = 6;
%
% for iroi = 1 : nrois
%     idx = a.dt.ecog.elec_roi{iroi};
%     [x{iroi}, y{iroi}, s{iroi}] =ecogFitPRFPix2Deg(subj, x0(idx), y0(idx), s0(idx));
%     ecc{iroi} = sqrt(x{iroi}.^2 + y{iroi}.^2);
% end
% param.elecEcc = ecc;
%
% dt = a.dt;
% dt.ecog.prf.x = x;
% dt.ecog.prf.y = y;
% dt.ecog.prf.s = s;
% dt.ecog.prf.ecc = ecc;
%
% save(fullfile(dataLoc, dt_fName), 'dt')

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

%% COMPUTE OFFSET RESPONSE INDEX:

% This is what I decided to do:
% (1) Use a parameter grid to find the best parameter set that I will use
%     as seed for the fine fit.
% (2) I do two kinds of fine fit here: First I fit the DN model to the
%     first 700ms of the time course; second I fit the DN model to the entire
%     time course.
% (3) Compare results from the two types of fine fit, and decide for each
%     electrode that parameters from which method I am going to rely on

% The reason I don't just use the model fit to the first 700ms of the data
% is that the model fit to the first 700ms of the time course do not
% necessarily predict the rest of the data. For the electrodes with
% relatively small offset response index, the model fit to the first 700ms
% of the data tend not to explain the rest of the data well.

%% INITIATE COMPUTING OFFSET RESPONSE INDEX

ORI = []; irfType = 'uniphasic'; trimIdx = 1 : 700; tTrim = t(trimIdx); stimTrim = stim(trimIdx);

%% COMPUTE OFFSEST RESPONSE INDEX - GRID FIT TO THE ENTIRE TIME COURSE

for iroi = whichRois
    nElec = length(ecc{iroi});
    toFit = bbts{iroi}';
    [ORI.seed{iroi}, ORI.seedR2{iroi}] = dn_gridFit(toFit, dn, stim, t, irfType);
end

%% COMPUTE OFFSET RESPONSE INDEX - GRID FIT TO THE FIRST 0.7 OF THE TIME COURSE

for iroi = whichRois
    nElec = length(ecc{iroi});
    toFit = bbts{iroi}(trimIdx, :)';
    [ORI.seedTrim{iroi}, ORI.seedTrimR2{iroi}] = dn_gridFit(toFit, dn, stimTrim, tTrim, irfType);
end

%% COMPUTE OFFSET RESPONSE INDEX - FIT THE ENTIRE TIME COURSE

for iroi = whichRois
    nElec = length(ecc{iroi});
    toFit = bbts{iroi}';
    seed  = [ORI.seed{iroi}, zeros(nElec, 1)];
    seedTrim = [ORI.seedTrim{iroi}, zeros(nElec, 1)];
    
    % COMPUTE MODEL FIT TO BOTH THE FIRST 700MS AND THE ENTIRE TIME COURSE
   % [param_fl{iroi}, ~, ~, exitflg_fl{iroi}] = dn_fineFit(toFit, stim, t, dn_ORI, seed, irfType);
    [param_pt{iroi}, ~, ~, exitflg_pt{iroi}] = dn_fineFit(toFit(:, trimIdx), stimTrim, tTrim, dn_ORI, seedTrim, irfType);
end

%% COMPARE BOTH WAYS OF COMPUTATION (PREDICTIONS AND R2)

select = {};
for iroi = whichRois
    nElec = length(ecc{iroi});
    data  = bbts{iroi};
    for k = 1 : nElec
        % COMPUTE MODEL PREDICTIONS ---------------------------------------
        thisprm_fl = dn_fillDNParams(param_fl{iroi}(k, :), irfType);
        thisprm_pt = dn_fillDNParams(param_pt{iroi}(k, :), irfType);
        % Use both sets of parameters to compare the prediction to the
        % first 700ms if the data:
        prd_fl_700{iroi}(k, :) = normMax(dn_DNmodel(thisprm_fl, stim(trimIdx), t(trimIdx)));
        prd_pt_700{iroi}(k, :) = normMax(dn_DNmodel(thisprm_pt, stim(trimIdx), t(trimIdx)));
        
        prd_fl_1204{iroi}(k, :) = normMax(dn_DNmodel(thisprm_fl, stim, t));
        prd_pt_1204{iroi}(k, :) = normMax(dn_DNmodel(thisprm_pt, stim, t));
        
        % COMPUTE R2 ------------------------------------------------------
        r2_fl{iroi}(k) = corr(data(trimIdx, k), prd_fl_1204{iroi}(k, trimIdx)').^2;
        r2_pt{iroi}(k) = corr(data(trimIdx, k), prd_pt_1204{iroi}(k, trimIdx)').^2;
        if r2_fl{iroi}(k) > r2_pt{iroi}(k), select{iroi}(k) = 1; else, select{iroi}(k) = 2; end
    end
end

%% MAKE ONE MINOR ADJUSTMENT BY HAND

% After observing figure (2), we need to make one adjustment by hand, the
% 11th electrode in V1.

% Should do this after the VISUALIZATION STEP (next), put it here only for
% convenince for saving figures.

% select{1}(11) = 1;

%% VISUALIZE TWO SETS OF MODEL PREDICTIONS TO DATA OF DIFFERENT LENGTHS

fg1 = figure (1); clf; fg2 = figure (2); clf;
for iroi = whichRois
    nElec = length(ecc{iroi});
    data  = bbts{iroi};
    for k = 1 : nElec
        data_700 = normMax(data(trimIdx, k)); data_1204 = normMax(data(:, k));
        figure (1)
        subplot(3, 15, (iroi - 1) * 15 + k),
        % PLOT DATA -------------------------------------------------------
        plot(t(trimIdx), data_700, 'k', 'linewidth', 2), hold on,
        % PLOT FULL MODEL PREDICTION --------------------------------------
        plot(t(trimIdx), prd_fl_700{iroi}(k, :), 'b-', 'linewidth', 1.5)
        % PLOT PARTIAL MODEL PREDICTION -----------------------------------
        plot(t(trimIdx), prd_pt_700{iroi}(k, :), 'r-', 'linewidth', 1.5)
        axis tight, box off, set(gca, 'ytick', [0, 1], 'xtick', [0.2, 0.7])
        if k == 1, legend('data', 'full', 'part', 'Location', 'best'), end
        
        figure (2)
        subplot(3, 15, (iroi - 1) * 15 + k),
        % COLOR THE WINNING MODEL -----------------------------------------
        if select{iroi}(k) == 1, patch([0, 1024, 1024, 0], [-0.1, -0.1, 1, 1], [0.85, 0.85, 1]); hold on
        else patch([0, 1024, 1024, 0], [-0.1, -0.1, 1, 1], [1, 0.85, 0.85]); hold on,
        end
        % PLOT DATA -------------------------------------------------------
        plot(t, data_1204, '-', 'linewidth', 2, 'color', 0.4 * ones(1, 3)),
        % PLOT FULL MODEL PREDICTION --------------------------------------
        plot(t, prd_fl_1204{iroi}(k, :), 'b-', 'linewidth', 1.8)
        % PLOT PARTIAL MODEL PREDICTION -----------------------------------
        plot(t, prd_pt_1204{iroi}(k, :), 'r-', 'linewidth', 1.8)
        axis tight, box off, set(gca, 'ytick', [0, 1], 'xtick', [0.2, 0.7]), xlim([0, 1]), ylim([-0.1, 1])
        
        if k == 1, legend('data', 'full', 'part', 'Location', 'best'), end
    end
end
fg1.Position = [1, 1000, 2500, 500]; fg2.Position = [1, 1000, 2500, 500];

%% SAVE FIGURES

printnice(2, 0, figureLoc, 'ana_chooseModelFitForORI');

%% COMPUTE THE OFFSET RESPONSE INDEX

ORI.rng = [701 : 1000]; ORI.index = {};

for iroi = whichRois
    nElec = length(ecc{iroi});
    data  = bbts{iroi};
    for k = 1 : nElec
        dt_max = normMax(data(:, k));
        dt_rng = dt_max(ORI.rng);
        
        if select{iroi}(k) == 1,
            prd_rng = prd_fl_1204{iroi}(k, ORI.rng);
            ORI.pred{iroi}(k, :) = prd_fl_1204{iroi}(k, :);
            
        else prd_rng = prd_pt_1204{iroi}(k, ORI.rng);
            ORI.pred{iroi}(k, :) = prd_pt_1204{iroi}(k, :);
        end
        ORI.index{iroi}(k) = mean(dt_rng - prd_rng');
    end
end

%% VISUALIZE THE FIT

fg3 = figure (3); clf
for iroi = whichRois
    for k = 1 : length(ecc{iroi})
        subplot(3, 15, (iroi - 1) * 15 + k)
        if elec_eccIdx{iroi}(k) == 1, patch([0, 1.24, 1.24, 0], [-0.1, -0.1, 1, 1], 0.6 * ones(1, 3)); hold on
        elseif elec_eccIdx{iroi}(k) == 2, patch([0, 1.24, 1.24, 0], [-0.1, -0.1, 1, 1], 0.85 * ones(1, 3)); hold on
        end
        
        plot(t, normMax(bbts{iroi}(:, k)), 'k-', 'linewidth', 2), hold on
        plot(t, ORI.pred{iroi}(k, :), 'r-', 'linewidth', 2), axis tight
        set(gca, 'xtick', [0.2, 0.7, 0.9]), ylim([-0.1, 1])
    end
end
fg3.Position = [1, 1000, 2500, 500];
printnice(3, 0, figureLoc, 'ana_visualizeModelFitForORI');

%% VISUALIZE THE OFFSET RESPONSE INDEX

figure (4), clf
for iroi = whichRois
    subplot(1, 3, iroi)
    idx = elec_eccIdx{iroi};
    for k = 1 : 3
    plot(k, mean(ORI.index{iroi}(idx == k)), 'ko', 'markerfacecolor', 'k', 'markersize', 10), hold on
    end
    %ylim([0, 0.15])
end

%% COMPUTE BIPHASIC FIT TO V1 - V3 ELECTRODES : GRID

%% COMPUTE BIPHASIC FIT TO V1 - V3 ELECTRODES : FINE FIT