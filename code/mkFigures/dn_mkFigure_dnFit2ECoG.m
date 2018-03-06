% dn_mkFigure_dnFit2ECoG

%% SAVE FIGURE KNOB

saveFigures = 1; % 0 if do not want to save figure

%% LOAD DATA AND PARAM FILES

dataLoc    = fullfile(dn_ECoG_RootPath, 'data');
dt_fName   = 'dn_preprocessedData.mat';
prm_fName  = 'dn_params.mat';

a = load(fullfile(dataLoc, dt_fName));
b = load(fullfile(dataLoc, prm_fName));

%% EXTRACT AND DERIVE PARAMETERS

ecogPrm = b.prm.ecog;

t    = ecogPrm.t;
stim = ecogPrm.stim;
dn   = ecogPrm.dn;
rois = ecogPrm.roiNm;

% DN MODEL PREDICTIONS TO THE ECOG DATA -----------------------------------
prd  = dn.prd;
dt   = dn.bs_bbts';
prm  = dn.param; % tau1, tau2, n, sigma, shift
metrics(1, :) = dn.summaryMetrics.t2pk; % time to peak
metrics(2, :) = dn.summaryMetrics.r_asymp;  % asymptotic response
r2   = dn.r2;

% DERIVED PARAMETERS ------------------------------------------------------
nBoots = size(dt, 1);
nRois  = length(rois);

for iroi = 1 : nRois
   nElec(iroi) = size(a.dt.ecog.bbts_roi{iroi}, 2); 
end

%% FIT SCALE TO THE ECOG DATA

% When we fit the DN model to the ECoG data, we normalize both the data and
% the model prediction to their maximum response. Here we fit the scale
% factor so that both the prediction and and data are in unit of percentage
% change.

irfType = 'uniphasic';

for iBoot = 1 : nBoots
    toFit = dt(iBoot, :);
    scale(iBoot) = dn_fitScale(toFit, prd(iBoot, :));
end

%% COMPUTE MEAN AND STD FOR MODEL RPEDICTIONS

for iBoot = 1 : nBoots
    prd(iBoot, :) = scale(iBoot) * normMax(prd(iBoot, :));
end

%% COMPUTE R2 RANGE FOR EACH ROI

r2_rng = [];

for iroi = 1 : nRois
    idx = (iroi - 1) * 100 + 1 : iroi * 100;
    r2_rng(iroi, :) = prctile(r2(idx), [25, 75]).*100;
end

%% COMPUTE MEAN AND STD FOR DATA AND MODEL PREDICTIONS

for iroi = 1 : nRois
    idx = (iroi - 1) * 100 + 1 : iroi * 100;
    % FOR DATA
    mdt(iroi, :) = mean(dt(idx, :));
    sdt(iroi, :) = std(dt(idx, :));
    % FOR MODEL PREDICTIONS
    mprd(iroi, :) = mean(prd(idx, :));
end

%% FIT DN MODEL TO THE MEAN DATA FOR EACH ROI

mseed = []; mprm = []; ext_mprm = zeros(nRois, 7);

% USE THE MEDIAN PARAMETERS AS SEED ---------------------------------------
for iroi = 1 : nRois
    idx  = (iroi - 1) * 100 + 1 : iroi * 100;
    mseed(iroi, :) = median(prm(idx, 1 : 5));
end

% FIT DN MODEL TO THE MEAN DATA -------------------------------------------
mprm = dn_fineFit(mdt, stim, t, ecogPrm.dn, mseed, irfType);

% COMPUTE MODEL PREDICTION WITHOUT SHIFT ----------------------------------
ext_mprm(:, [1, 3 : 5]) = mprm(:, 1 : 4);
ext_mprm(:, 7)          = 1;
for iroi = 1 : nRois
    mprm_prd(iroi, :) = normMax(dn_DNmodel(ext_mprm(iroi, :), stim, t));
end

% COMPUTE MODEL PREDICTION WITH SHIFT -------------------------------------
sft_mprm(:, [1, 3 : 6]) = mprm(:, 1 : 5);
sft_mprm(:, 7)          = 1;
for iroi = 1 : nRois
    msft_prd(iroi, :) = normMax(dn_DNmodel(sft_mprm(iroi, :), stim, t));
end

%% COMPUTE MEAN AND STD OF THE SUMMARY METRICS

m_metrics = []; s_metrics = [];

for iroi = 1: nRois
    idx = (iroi - 1) * 100 + 1 : iroi * 100;
    m_metrics(iroi, :)    = median(metrics(:, idx), 2);
    s_metrics(iroi, :, :) = prctile(metrics(:, idx), [25, 75], 2);
end

%% PLOT ECOG BROADBAND DATA AND MODEL PREDICTIONS (PANEL A)

fg1 = figure (1); clf
for iroi = 1 : nRois
    subplot(2, 3, iroi)
    % plot stimulus
    patch([0, 0.5, 0.5, 0], [0, 0, 7.6, 7.6], 0.95 * ones(1, 3)), hold on
    % plot data and model predictions
    shadedErrorBar(t-0.2, mdt(iroi, :), sdt(iroi, :), 'b'),
    plot(t-0.2, mprd(iroi, :), 'r-', 'linewidth', 2)
    axis tight, ylim([-0.2, 7.6]), box off,
    set(gca, 'xaxislocation', 'origin', 'fontsize', 16, 'xtick', [0, 0.5]), title(rois{iroi})
    % number of electrodes per ROI
    text(0.6, 5, sprintf('n = %d', nElec(iroi)), 'fontsize', 16), 
    % r2 ramge
    text(0.6, 4, sprintf('r^2 : [%.1f, %.1f]', r2_rng(iroi, 1), r2_rng(iroi, 2)), 'fontsize', 16)
end
subplot(2, 3, 1), ylabel('fractional signal change'), xlabel('time (s)')

fg1.Position = [1, 500, 2000, 500];

%% COMPARE RESPONSE SHAPES ACROSS ROIS (PANEL B)

% FIRST, SANITY CHECK, TO CONFIRM THE FIT TO THE MEAN DATA TIME COURSE IS
% GOOD
figure (100), clf
for iroi = 1 : nRois
    subplot(2, 3, iroi)
    plot(t-0.2, normMax(mdt(iroi, :)), 'k-', 'linewidth', 3), hold on,
    plot(t - 0.2, msft_prd(iroi, :), 'r-', 'linewidth', 2)
    plot(t - 0.2, mprm_prd(iroi, :), 'b-', 'linewidth', 2), % predicted ts without shift
    set(gca, 'xaxislocation', 'origin'), box off, axis tight, ylim([-0.1, 1]), set(gca, 'fontsize', 14)
    title(rois{iroi}), set(gca, 'xtick', [0, 0.5])
end
subplot(2, 3, 1), legend('data', 'prd. with shift', 'prd without shift'),

% COMPARE RESPONSE SHAPES ACROSS ROIS -------------------------------------
fg2 = figure (2); clf, colorMatrix = copper(6);
for k = 1 : 2
    subplot(1, 2, k),
    patch([0, .5, .5, 0], [0, 0, 1, 1], 0.9 * ones(1, 3)), hold on
    for iroi = 1 : nRois
        plot(t - 0.2, mprm_prd(iroi, :), 'color', colorMatrix(iroi, :), 'linewidth', 2),
        box off, axis tight, set(gca, 'ytick', [0, 0.5, 1], 'xtick', [0, 0.5])
    end
    legend(['stim', rois]), set(gca, 'fontsize', 14)
end
subplot(1, 2, 2), xlim([0, 0.3]), ylim([0.3, 1]), axis square

fg2.Position = [1, 250, 1000, 250];

%% PLOT SUMARY METRICS (PANEL C)

title_txt = {'T_{peak}', 'R_{asymp}'};

fg3 = figure (3); clf
for k = 1 : 2
    subplot(1, 2, k)
    % PLOT THE MEDIAN OF THE ESTIMATED SUMMARY METRICS
    plot(m_metrics(:, k), 'ko', 'markerfacecolor', 'k', 'markersize', 6), hold on
    % PLOT THE ERROR BARS OF THE ESTIMATED SUMMARY METRICS
    for iroi = 1 : nRois
        plot([iroi, iroi], squeeze(s_metrics(iroi, k, :)), 'k-')
    end
    % OTHER FIGURE FEATURES
    xlim([0.5, nRois + 0.5]), axis square, box off,
    set(gca, 'xtick', 1 : nRois, 'xticklabel', rois, 'fontsize', 16), title(title_txt{k})
end
subplot(1, 2, 1), ylim([0, 180]), subplot(1, 2, 2), ylim([0, 0.2])

fg3.Position = [1, 500, 1000, 500];
%% SAVE FIGURES

if saveFigures
    figLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
    % PANEL A
    fg1Nm = 'fg_dnFit2ECoG_A';
    printnice(1, 0, figLoc, fg1Nm);
    % PANEL B
    fg2Nm = 'fg_dnFit2ECoG_B';
    printnice(2, 0, figLoc, fg2Nm);
    % PANEL C
    fg3Nm = 'fg_dnFit2ECoG_C';
    printnice(3, 0, figLoc, fg3Nm);
end

%% 

figure 
for k = 1 : 4
   subplot(1, 4, k)
   plot(mprm(:, k), 'ko')
end
