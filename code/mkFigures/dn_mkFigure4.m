% dn_mkFigure4

%% LOAD DATA AND PARAM FILES

dataLoc    = fullfile(dn_ECoG_RootPath, 'data');
dt_fName   = 'dn_preprocessedData.mat';
prm_fName  = 'dn_params.mat';

a = load(fullfile(dataLoc, dt_fName));
b = load(fullfile(dataLoc, prm_fName));

%% EXTRACT PARAMETERS

dn   = b.prm.ecog.dn;

metrics(1, :) = dn.summaryMetrics.t2pk; % time to peak
metrics(2, :) = dn.summaryMetrics.r_asymp;  % asymptotic response

prm  = dn.param; % tau1, tau2, n, sigma, shift

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

%% PLOT SUMMARY METRICS

title_txt = {'T_{peak}', 'R_{asymp}'};

fg3 = figure (1); clf
for k = 1 : 2
    subplot(4, 4, k)
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
subplot(4, 4, 1), ylim([100, 160]), subplot(4, 4, 2), ylim([0, 0.14])

%% PLOT MODEL PARAMETERS

title_txt = {'tau1', 'tau2', 'n', 'sigma'};

for k = 1 : 4
    subplot(4, 4, k+4), cla
    plot(m_prm(:, k), 'ko', 'markerfacecolor', 'k'), hold on
    for iroi = 1 : nRois
        plot([iroi, iroi], squeeze(s_prm(iroi, :, k)), 'k-')
    end
    xlim([0.5, nRois + 0.5]), axis square, box off,
    set(gca, 'xtick', 1 : nRois, 'xticklabel', rois, 'fontsize', 16),title(title_txt{k})
end
subplot(4, 4, 5), ylim([0, 0.45])
subplot(4, 4, 6), ylim([0, 0.2])
subplot(4, 4, 7), ylim([0, 5])
subplot(4, 4, 8), ylim([0, 0.2])
