% dn_ecogSubj


saveFigure = 1;

%% SUBJECT INFORMATION

% Five runs

% sub-19 107: V3v?
% sub-19 108: V2v?
% sub-19 109: V1?
% sub-19 115: V2d?
% sub-19 120: V3AB?
% sub-19 121: V3AB?

whichElec = [107, 108, 109, 115, 120, 121];


%% USEFUL FUNCTIONS

whiten         = @(x) (x - mean(x(:)))./ diff(prctile(x, [.25 .75]));
whiten_hilbert = @(x) abs(hilbert(whiten(x)));

%% KNOBS

visualizeStim = 1;

%% STIMULUS AND DATA LOCATIONS

data = [];

dataLoc = fullfile('/Volumes', 'server', 'Projects', 'ECoG', 'soc');
subjPth = fullfile(dataLoc, 'sub-19');
stimPth = fullfile(dataLoc, 'stimuli');
dataPth = fullfile(subjPth, 'ses-01', 'ieeg');
stimNm  = 'task-soc_stimuli.mat';

dtFiles = dir(fullfile(dataPth, 'sub-19*.mat'));
eventFiles = dir(fullfile(dataPth, '*events*'));

%% LOAD MODEL PARAMETERS

prmLoc = fullfile(dn_ECoG_RootPath, 'data');
fNm    = 'dn_params.mat';

c = load(fullfile(prmLoc, fNm));

ecogPrm = c.prm.ecog.dn;

%% DERIVED PARAMETERS
raw   = [];

nRuns = 5;
nElec = length(whichElec);

%% LOAD DATA, STIMULUS AND EVENT FILES

% LOAD STIMULUS FILE
a = load(fullfile(stimPth, stimNm));

% LOAD DATA FILES
for k = 1 : nRuns
    b = load(fullfile(dataPth, dtFiles(k).name));
    data{k} = b.data(:, whichElec);
end

srate = b.srate;
dt    = 1/srate;

% LOAD EVENT FILES
for k = 1 : nRuns
    c = fopen(fullfile(dataPth, eventFiles(k).name));
    C = textscan(c, '%f %f %d %s %d', 'HeaderLines', 1);
    raw.trialStart{k} = C{1}(1 : 2 : end);
    raw.stimType{k}   = C{3}(1 : 2 : end);
end

%% VISUALIZE STIMULI

% SELECT FULLFIELD/FULLFIELD NOISE STIMULI --------------------------------
fullfield_idx      = [10, 29, 39 : 86];
fullfield_noiseIdx = [10, 29, 55 : 68, 74 : 78, 83 : 86];

if visualizeStim
    figure (1), clf, colormap gray
    for k = 1 : 86
        subplot(9, 10, k)
        imagesc(squeeze(a.stimuli(:, :, k))), axis off
        if ismember(k, fullfield_noiseIdx), title(k); end
    end
end

%% CREATE TRIAL TIMING INDEX

for irun = 1 : nRuns
    nTrials    = length(raw.trialStart{irun});
    thisTrSt   = round(raw.trialStart{irun}, 3);
    t          = round(dt : dt : dt * size(data{irun}, 1), 3);
    % make trial idx
    for iTrial = 1 : nTrials
        st = find(t == thisTrSt(iTrial));
        tIdx{irun}(iTrial, :) = st-199 : st + 1000;
    end
end

%% CREATE STIMULUS INDEX

for irun = 1 : nRuns
    relevantStim{irun} = ismember(raw.stimType{irun}, fullfield_noiseIdx);
end

%% EXTRACT BROADBAND

%% BAND-PASS FILTER THE TIME COURSE FOR EACH RUN

bands = [70 80; 80 90; 90 100; 100 110; 130 140; 140 150; 150 160; 160 170; 190 200; 200 210];

% BANDPASS FILTER THE ORIGINAL DATA ---------------------------------------
for irun = 1 : nRuns
    bp{irun} = bandpassFilter(data{irun}, srate, bands); % time course x n electrodes x 10 bands
end

raw.bp = bp;

%% EXTRACT BROADBAND USING THE BANDPASS FILTERED DATA

bb_elec = {}; bb_rs = {};

for irun = 1 : nRuns
    for iElec = 1 : nElec
        bp_elec = squeeze(bp{irun}(:, iElec, :));
        bb_elec = sum(whiten_hilbert(bp_elec).^2, 2);
        
        % RESHAPE THE EXTRACTED BROADBAND DATA
        for iTrial = 1 : nTrials
            bb_rs{irun}(iTrial, :, iElec) = bb_elec(tIdx{irun}(iTrial, :));
            %reshape(bb_elec, size(tIdx{irun}));
        end
    end
end

%% AVERAGE OVER TRIALS

raw.bb = [];

for irun = 1 : nRuns
    raw.bb(irun, :, :) = squeeze(mean(bb_rs{irun}(relevantStim{irun}, :, :)));
end

raw.mbb = squeeze(mean(raw.bb));

%% VISUALIZE

order = [3, 2, 4, 1, 5, 6];

for iElec = 1 : nElec
    baseline = mean(raw.mbb(1 : 200, iElec))
    raw.mbb(:, iElec) = raw.mbb(:, iElec) ./baseline;
    raw.mbb(:, iElec) = raw.mbb(:, iElec) - 1;
    ts_max(iElec) = max(raw.mbb(:, iElec));
end

trial_t = -.2 + dt : dt : 1;

figure (2), clf,
subplot(2, 2, 1), set(gca, 'colororder', copper(6)), hold on
patch([0, 0.5, 0.5, 0], [0, 0, 1, 1], 0.9 * ones(1, 3))
plot(trial_t, raw.mbb(:, order ), 'linewidth', 2),
legend({'stim', 'V3v?', 'V2v?', 'V1?', 'V2d?', 'V3AB?', 'V3AB?'});
set(gca, 'fontsize', 14, 'xaxislocation', 'origin')

subplot(2, 2, 2), set(gca, 'colororder', copper(6)), hold on
patch([0, 0.5, 0.5, 0], [0, 0, 1, 1], 0.9 * ones(1, 3))
plot(trial_t, normMax(raw.mbb(:, order)), 'linewidth', 2),
%legend({'stim', 'V3v?', 'V2v?', 'V1?', 'V2d?', 'V3AB?', 'V3AB?'});
legend({'stim', 'V1', 'V2v', 'V2d', 'V3v', 'v3AB', 'V3AB'})
set(gca, 'fontsize', 14, 'xaxislocation', 'origin')
xlim([0, 0.3]), ylim([0.5, 1])

%% GROUP THE BROADBAND DATA

group     = {3, [2,4], [1, 5, 6]}; % V1, V2, and V3 +

for k = 1 : length(group)
    raw.tofit(k, :) = normMax(mean(normMax(raw.mbb(:, group{k})), 2));
    %raw.tofitNoNorm(k, :) = mean(raw.mbb(:, group{k}), 2);
    tmp = mean(raw.mbb(:, group{k}));
    group_max(k) = mean(ts_max(group{k}));
end

fg2 = figure (2);
for k = 1 : 2
    subplot(2, 2, k + 2), cla, set(gca, 'colororder', copper(3)), hold on
    patch([0, .5, .5, 0], [0, 0, 1, 1], 0.9 * ones(1, 3))
    plot(trial_t, normMax(raw.tofit'), 'linewidth', 3), set(gca, 'xaxislocation', 'origin')
    set(gca, 'fontsize', 14), set(gca, 'xtick', [0, 0.5], 'ytick', [0, 1]), xlabel('time (s)')
    legend('stimulus', 'V1', 'V2', 'V3 and anterior'), axis tight
end
subplot(2, 2, 4), xlim([0.05, 0.21]), ylim([0.5, 1])

fg2.Position = [500, 800, 800, 500];

%% FIT THE DN MODEL TO THE DATA

%  dn_fineFit(data, stim, t, param, seed, irfType)

seed = repmat([0.07, 0.07, 3, 0.05, 0], [3, 1]);
t    = [1 : length(trial_t)]./1000;
stim = zeros(1, length(t));
stim(200 : 700) = 1;

irfType = 'uniphasic';
raw.param  = dn_fineFit(raw.tofit, stim, t, ecogPrm, seed, irfType);
derivedPrm = dn_computeDerivedParams(raw.param, irfType);

fg3 = figure (3);
subplot(1, 2, 1)
plot(derivedPrm.t2pk, 'ro', 'markersize', 7, 'markerfacecolor', 'r'), xlim([0.5, 3.5]), ylim([100, 120]), axis square, box off
set(gca, 'xtick', 1 : 3, 'xticklabel', {'V1', 'V2', 'V3/anterior'}), set(gca, 'fontsize', 14), title('T_{peak} (ms)')

subplot(1, 2, 2)
plot(derivedPrm.r_asymp, 'ro', 'markersize', 7, 'markerfacecolor', 'r'), xlim([0.5, 3.5]), ylim([0.14, 0.18]), axis square, box off
set(gca, 'xtick', 1 : 3, 'xticklabel', {'V1', 'V2', 'V3/anterior'}), set(gca, 'fontsize', 14), title('R_{asymp}')

fg3.Position = [250, 500, 500, 250];

%% PLOT THE MODEL FIT

% Generate model predictions

for k = 1 : 3
   ext_prm = [raw.param(k, 1), 0, raw.param(k, 2 : 5), 1];
   raw.prd(k, :) = normMax(dn_DNmodel(ext_prm, stim, t));
end

fg4 = figure (4); clf, title_txt = {'V1', 'V2', 'V3/anterior'};

for k = 1 : 3
    subplot(1, 3, k)
    patch([0, .5, .5, 0], [0, 0, 1, 1], 0.9 * ones(1, 3)),  hold on
    plot(trial_t, raw.tofit(k, :), 'k', 'linewidth', 2), axis tight
    plot(trial_t, raw.prd(k, :), 'r-', 'linewidth', 2), set(gca, 'xaxislocation', 'origin'), set(gca, 'xtick', [0, 0.5], 'ytick', [0, 1])
    set(gca, 'fontsize', 14)
    title(title_txt{k})
end

fg4.Position = [250, 2000, 2000, 250];

% PLOT MODEL PREDICTIONS ALL TOGETHER
fg5 = figure (5), clf

for k = 1 : 3
    ext_prm = [raw.param(k, 1), 0, raw.param(k, 2 : 4),0, 1];
    prd(k, :)     = normMax(dn_DNmodel(ext_prm, stim, t));
end

for k = 1 : 2
    subplot(1, 2, k), set(gca, 'colororder', copper(3)), hold on
    patch([0, .5, .5, 0], [0, 0, 1, 1], 0.9 * ones(1, 3))
    plot(trial_t, prd, '-', 'linewidth', 3), set(gca, 'xaxislocation', 'origin'), set(gca, 'xtick', [0, 0.5], 'ytick', [0, 1])
    set(gca, 'xaxislocation', 'origin')
    set(gca, 'fontsize', 14), set(gca, 'xtick', [0, 0.5], 'ytick', [0, 1]), xlabel('time (s)')
    legend('stimulus', 'V1', 'V2', 'V3 and anterior'), axis tight
end
subplot(1, 2, 2), xlim([0.05, 0.21]), ylim([0.5, 1])

%% SAVE FIGURE

if saveFigure, 
    figLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
    fg5Nm = 'fg_additionalData_A';
    printnice(5, 0, figLoc, fg5Nm);
    
    fg3Nm = 'fg_additionalData_B';
    printnice(3, 0, figLoc, fg3Nm);
    
    fg4Nm = 'fg_additionalData_C';
    printnice(4, 0, figLoc, fg4Nm);
end


