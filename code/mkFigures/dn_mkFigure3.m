% dn_mkFigure3

%% LOAD DATA AND PARAM FILES

dataLoc    = fullfile(dn_ECoG_RootPath, 'data');
dt_fName   = 'dn_preprocessedData.mat';
prm_fName  = 'dn_params.mat';

a = load(fullfile(dataLoc, dt_fName));
b = load(fullfile(dataLoc, prm_fName));

%% EXTRACT AND DERIVE PARAMETERS

irfType = 'uniphasic';

ecogPrm = b.prm.ecog;

t    = ecogPrm.t;
stim = ecogPrm.stim;
dn   = ecogPrm.dn;
rois = ecogPrm.roiNm;
orig_dt = a.dt.ecog.bbts_roi;

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

% % When we fit the DN model to the ECoG data, we normalize both the data and
% % the model prediction to their maximum response. Here we fit the scale
% % factor so that both the prediction and and data are in unit of percentage
% % change.

% compute the max level of the original data in each ROI:
for iroi = 1 : nRois
    t_max           = max(orig_dt{iroi});
    mt_max(iroi, :) = mean(t_max);
end

scale = repmat(mt_max, [1, 100])';
scale = reshape(scale, [nRois * 100, 1]);

%% COMPUTE MEAN AND STD FOR MODEL RPEDICTIONS

for iBoot = 1 : nBoots
    prd(iBoot, :) = scale(iBoot) * normMax(prd(iBoot, :));
    dt(iBoot, :)  = scale(iBoot) * normMax(dt(iBoot, :));
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

%% COMPUTE THE MEDIAN AND INTERQUARTILE RANGE FOR THE MODEL PARAMETERS

m_prm = []; s_prm = [];

for iroi = 1 : nRois
    idx = (iroi - 1) * 100 + 1 : iroi * 100;
    m_prm(iroi, :) = median(prm(idx, :));
    s_prm(iroi, :, :) = prctile(prm(idx, :), [25, 75]);
end

%% PLOT ECOG BROADBAND DATA AND MODEL PREDICTIONS (PANEL A)

fg1 = figure (1); clf
for iroi = 1 : nRois
    subplot(2, 4, iroi)
    % plot stimulus
    patch([0, 0.5, 0.5, 0], [0, 0, 7.6, 7.6], 0.95 * ones(1, 3)), hold on
    % plot data and model predictions
    shadedErrorBar(t-0.2, mdt(iroi, :), sdt(iroi, :), 'k'),
    plot(t-0.2, mprd(iroi, :), 'r-', 'linewidth', 2)
    axis tight, ylim([-0.2, 7.6]), box off,
    set(gca, 'xaxislocation', 'origin', 'fontsize', 16, ...
        'xtick', [0, 0.5], 'YLim', [-0.1 7]), title(rois{iroi})
    % number of electrodes per ROI
    text(0.6, 5, sprintf('n = %d', nElec(iroi)), 'fontsize', 16),
    % r2 ramge
    text(0.6, 4, sprintf('r^2 : [%.1f, %.1f]', r2_rng(iroi, 1), r2_rng(iroi, 2)), 'fontsize', 16)
end
subplot(2, 4, 1), ylabel('fractional signal change'), xlabel('time (s)')

fg1.Position = [1, 300, 2000, 300];

%% COMPARE RESPONSE SHAPE ACROSS ROIS

% COMPARE RESPONSE SHAPES ACROSS ROIS -------------------------------------
colorMatrix = copper(nRois);
for k = 1 : 2
    subplot(2, 4, 4+k),
    patch([0, .5, .5, 0], [0, 0, 1, 1], 0.9 * ones(1, 3)), hold on
    for iroi = 1 : nRois
        plot(t - 0.2, mprm_prd(iroi, :), 'color', colorMatrix(iroi, :), 'linewidth', 2),
        box off, axis tight, set(gca, 'ytick', [0, 0.5, 1], 'xtick', [0, 0.5])
    end
    legend(['stim', rois]), set(gca, 'fontsize', 14)
end
subplot(2, 4, 5),  xlim([0, 0.3]), ylim([0.3, 1]), axis square

%% PLOT CROSS-VALIDATED RESULTS (across trials)

trial_xr2 = b.prm.ecog.dn.xr2.trial_r2;

subplot(2, 4, 7), cla
for iroi = 1 : 4
   plot(iroi, median(trial_xr2{iroi}), 'ro', 'markerfacecolor', 'r', 'markersize', 15), hold on 
   plot(iroi, trial_xr2{iroi}, 'ko', 'markerfacecolor', 'k', 'markeredgecolor', 'w'), 
end
box off, ylim([0, 1]), set(gca, 'fontsize', 16), xlim([0.5, 4.5]), set(gca, 'xtick', 1 : 4), 
set(gca, 'xticklabel', {'V1', 'V2', 'V3', 'anterior'}), title('x-validated R2 over trials')

%% PLOT CROSS-VALIDATED RESULTS (across electrodes)

elec_xr2 = b.prm.ecog.dn.xr2.elec_r2;

subplot(2, 4, 8), cla
for iroi = 1 : 4
   plot(iroi, median(elec_xr2{iroi}), 'ro', 'markerfacecolor', 'r', 'markersize', 15), hold on 
   plot(iroi, elec_xr2{iroi}, 'ko', 'markerfacecolor', 'k', 'markeredgecolor', 'w'), 
end
box off, ylim([0, 1]), set(gca, 'fontsize', 16), xlim([0.5, 4.5]), set(gca, 'xtick', 1 : 4), 
set(gca, 'xticklabel', {'V1', 'V2', 'V3', 'anterior'}), title('x-validated R2 over electrodes')