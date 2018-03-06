function [] = dn_spectrogram_spectra(stimulusGrouping, fgNum)
% function [] = dn_spectrogram()
%
% DESCRIPTION -------------------------------------------------------------
% This file makes a spectrogram for stimulus triggered response (0 - 500ms) after stimulus onsets) for every
% electrode. 
%
% INPUTS ------------------------------------------------------------------
% stimulusGrouping: 
% (1) 'all' means the computed spectrogram is averaged across all stimuli
% (2) 'noise' means the spectorgram is averaged across all noise stimuli ([1-3])
% (3) 'grating' means spectrogram averaged across all grating stimuli
% (4) 'k' means the spectrogram averaged across the k'th stimulus type, k
%      is an interger between 1 and 7.

% fgNum: figure Number

% DEPENDENCIES ------------------------------------------------------------
% knkutils: printnice.m
% Chronux functions for making spectrograms

%% EXAMPLE

stimulusGrouping = 'noise';
fgNum = 100;

%% SAVE FIGURE KNOB

saveFigure = 1;

%% IMPORT AND EXTRACT DATA

disp('dn_spectrogram: extracting data...')

fName = 'dn_rawData.mat';
dtLoc = fullfile(dn_ECoG_RootPath, 'data');
a     = load(fullfile(dtLoc, fName));
raw   = a.raw;

ts    = raw.ts;
idx   = raw.idx;
srate = 1000;

%% DERIVED PARAMETERS

nElec   = size(ts, 2); % number of electrodes
stimNms = raw.stimNames;

%% RE-SHAPE THE DATA INTO N_TRIALS X X_TIMEPOINTS_PER_TRIAL

tsRs = [];

% RE-SHAPE THE TIME COURSE
for k = 1 : nElec, tsElec = ts(:, k); tsRs(:, :, k) = tsElec(idx); end

%% AVERAGE ACROSS STIMULUS TYPES

switch stimulusGrouping
    case 'all'
        % do nothing, tsRs maintains its original shape
    case 'noise'
        tsRs = tsRs(stimNms < 4, :, :);
    case 'grating'
        tsRs = tsRs(stimNms > 3, :, :);
    case {1, 2, 3, 4, 5, 6, 7}
        tsRs = tsRs(stimNms == stimulusGrouping, :, :);
    otherwise
        error('The selected stimulusGrouping does not exist.')
end

%% DEFINE TIME-FREQUENCY MULTITAPER SETTINGS

% same setting as in Hermes et al. 2015, except that we average data across
% 3 stimulus types
params = [];

movingwin       = [.200 .05];
params.pad      = -1;
params.tapers   = [3 5];
params.fpass    = [0 200];
params.Fs       = srate;
params.trialave = 1;

baseline_tRange = [751 : 1200];
stimTrig_tRange = [251 : 700];

%% CAUCULATE THE SPECTROGRAM

disp('dn_spectrogram: computing the spectrogram...')

S_b = []; S_s = []; S = [];

for k1 = 1 : nElec
    for k = 1 : size(tsRs, 1)
        % BASELINE ------------------------------------------------------------
        [S_b(k1, k, :, :), t_b, f_b] = mtspecgramc(squeeze(tsRs(k, baseline_tRange, k1)), movingwin, params);
        
        % STIMULUS-TRIGGERED --------------------------------------------------
        [S_s(k1, k, :, :), t_s, f_s] = mtspecgramc(squeeze(tsRs(k, stimTrig_tRange, k1)), movingwin, params);
    end
end
% NORMALIZE WITH BASELINE ---------------------------------------------
S = squeeze(mean(S_s./S_b, 2));

%% VISUALIZE THE SPECTROGRAM

fg = figure (fgNum); clf, fg.Position = [1, 2000, 2000, 2000];

for k = 1 : nElec
    subplot(8, 10, k)
    imagesc(t_s, f_s, log10(squeeze(S(k, :, :)))), hold on, 
    axis xy
    set(gca,'XTick',[0 .4]), title(raw.goodChannels(k)), set(gca, 'fontsize', 12)
end

%% COMPUATE SPECTRA

f = [0 : length(baseline_tRange)-1]./(length(baseline_tRange)/1000);

spectra_stim = []; mspectra_stim = []; spectra_base = []; mspectra_base = [];

n = size(tsRs,2)/2;

for k = 1 : size(tsRs, 1)
   for iElec = 1 : nElec
      thisTS_stim = squeeze(tsRs(k, stimTrig_tRange, iElec));
      spectra_stim(k, :, iElec) = abs(fft(thisTS_stim))/n;
      
      thisTS_base = squeeze(tsRs(k, baseline_tRange, iElec));
      spectra_base(k, :, iElec) = abs(fft(thisTS_base))/n;
   end
end
% AVERAGE THE FREQUENCY COMPONENTS ACROSS STIMULI -------------------------
mspectra_stim = squeeze(mean(spectra_stim));
mspectra_base = squeeze(mean(spectra_base));

%% VISUALIZE SPECTRA

fg = figure (fgNum + 1), clf,  fg.Position = [1, 2000, 2000, 2000];
for iElec = 1 : nElec
    subplot(8, 10, iElec)
    loglog(f, mspectra_stim(:, iElec), 'k-', 'linewidth', 3), hold on, 
    loglog(f, mspectra_base(:, iElec), '-', 'linewidth', 3, 'color', 0.6 * ones(1, 3)),
    set(gca, 'xtick', [50, 100, 200]), ylim(10.^[-1 1]); title(raw.goodChannels(iElec)), set(gca, 'fontsize', 12), 
    xlim([30, 200]), box off
end

%% SAVE THE SPECTROGRAM IMAGES

if saveFigure,
    saveLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
    if isstr(stimulusGrouping)
        saveNm1  = sprintf('pre_spectrogram_indiElectrodes_%s', stimulusGrouping);
        saveNm2  = sprintf('pre_spectra_indiElectrodes_%s', stimulusGrouping);
    elseif isnumeric()
        saveNm1  = sprintf('pre_spectrogram_indiElectrodes_%d', stimulusGrouping);
        saveNm1  = sprintf('pre_spectra_indiElectrodes_%d', stimulusGrouping);
    end
    printnice(fgNum, 0, saveLoc, saveNm1);
    printnice(fgNum + 1, 0, saveLoc, saveNm2);
end

end
