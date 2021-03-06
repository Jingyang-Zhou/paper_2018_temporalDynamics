% dn_extract broadband

% WHAT IS STORED IN THE RAW DATA ------------------------------------------

% The ECoG experiment consists of 
% (1) 7 stimulus types, each type with 30 repeats, so a total of 210 trials
% (2) each trials consists of 1204 sampling point (time course, 1000Hz)
% (3) a total of 71 good channels and their correpsonding labels
% (4) an extracted broadband matrix, which we will compute here

%% LOAD, MAKE AND SAVE RAW DATA

raw = [];

[raw.ts, raw.stimNames, raw.imNames, raw.goodLabels, raw.goodChannels,raw.idx] = dn_getRawData;

% save raw data
% fName = 'dn_rawData.mat'; fLoc  = fullfile(dn_ECoG_RootPath, 'data', fName);
% save(fLoc, 'raw')

%% PRE-DEFINED VARIABLES AND FUNCTIONS

% USEFUL FUNCTIONS --------------------------------------------------------
whiten   = @(x) (x - mean(x(:)))./ diff(prctile(x, [.25 .75]));
normMax  = @(x) x./max(x);
normBase = @(x) x-mean(x(1 : 200));

% PREDEFINED VARIABLES ----------------------------------------------------
srate    = 1000;
nElec    = size(raw.ts, 2); % number of electrodes
exp_elec = 1;               % example electrode
exp_ts   = raw.ts(:, exp_elec);

% HILBERT TRANSFORMS ------------------------------------------------------
stand_hilbert = @(x) abs(hilbert(x));

%% ABOUT BROADBAND SIGNAL

% Broadband signal refers to a uniform increase (or a shift) of all the
% frequency bands within some range due to, say, a stimulus triggered
% neuronal firing rate increase.

% HOW TO EXTRACT BROADBAND? -----------------------------------------------
% The point of extracting broadband here is to look at the time course of,
% say, the increase in neuronal firing rate due to stimulus.

% This code is attempting to answer the following questions, some harder
% than the others

% 1. The bandpass filter
%    1.1 What range of frequencies should we bandpass filters?
%    1.2 How many bands should we use?
%    1.3 Are the answers to the above question universal?

% 2. The difference between averaging the bandpassed signals over trials
%    before Hilbert transform, or Hilbert transform before averaging.

%% QUESTION 1: HOW MANY BANDS SHOULD WE USE AND WHAT RANGE OF BANDS

bands{1} = {[60, 200], 140};
bands{2} = {[60, 200], 70};
bands{3} = {[60, 200], 20};
bands{4} = {[60, 200], 10};

% WHAT HAPPENS WHEN WE USE A DIFFERENT RANGE OF BANDS ---------------------

% The reason why we shouldn't include lower band is easier: the uniform
% shift in frequency amplitudes due to stimulus-triggered response does not
% apply to lower bands. The reason why we don't use only higher frequency
% bands is not known except for the obvious that we should use all the
% information that is valid.

bands{5} = {[10, 30], 20};
bands{6} = {[100, 200], 20};

% EXTRACT AND COMPARE BROADBAND USING THE ABOVE BANDWIDTH ----------------
bp = {}; bb = {}; m_bb = [];

for k = 1 : length(bands)
    nBands     = (bands{k}{1}(2) - bands{k}{1}(1))/bands{k}{2};
    bp{k}      = squeeze(bandpassFilter(exp_ts, srate, bands{k}));
    % take Hilbert transform 
    bb{k}      = geomean(stand_hilbert(bp{k}).^2, 2);
    m_bb(k, :) = mean(bb{k}(raw.idx));
end

%% VISUALIZE THE ANSWER TO QUESTION 1

figure (1), clf

for k = 1 : length(bands)
    subplot(length(bands), 1, k)
    plot(m_bb(k, :)), axis tight
end

% 
% 
% %% (1) BANDPASS ([60, 200]) FILTER THE SIGNAL THEN LOOK AT THE TIME COURSE
% 
% % So it seems like the first order thing to do is to look at the time
% % course of the band-passed signal:
% 
% bands  = {[40, 60], 20};
% nBands = (bands{1}(2) - bands{1}(1))/bands{2};
% 
% % band-pass filter the raw time series
% bp       = bandpassFilter(raw.ts, srate, bands);
% exp_bp   = bp(:, exp_elec); % look at an example electrode
% exp_bp   = exp_bp(raw.idx);
% m_exp_bp = mean(exp_bp);
% 
% % hilbert-transform (intuition: compute spike rate rather than looking at the raw signals)
% exp_bb = abs(hilbert(exp_bp));
% m_exp_bb = mean(exp_bb, 1);
% 
% % visualize 
% figure (1), clf
% %subplot(2, 2, 1), plot(m_exp_bp, 'k-'), hold on, plot(m_exp_bb, 'r-', 'linewidth', 2), axis tight
% 
% %% BANDPASS FILTER THE SIGNALS
% 
% % Variables to change:
% bands        = {[60, 200], 20}; % the bands are 60-80, 80-100, 100-120, 120-140,140-160,160-180,180-200
% nbands       = 7;
% 
% example_elec = 1;
% 
% nElec = size(raw.ts, 2);
% bp    = [];
% 
% % bandpass-filter the signal
% bp = bandpassFilter(raw.ts, srate, bands);
% 
% % make an example band-passed signal
% example_bp   = [];
% 
% for iband = 1 : nbands
%     tmp = squeeze(bp(:, example_elec, iband));
%     example_bp(:, :, iband) = tmp(raw.idx);
% end
% 
% % VISUALIZE BAND-PASS FILTER ---------------------------------------------
% 
% m_examp_bp = squeeze(mean(example_bp));
% 
% figure (1), clf
% for k = 1 : nbands
%     subplot(8, 1, k)
%     plot(m_examp_bp(:, k), 'k-'), hold on,
%     plot(abs(hilbert(m_examp_bp(:, k))), 'r-', 'linewidth', 2)
%     %plot(geomean(abs(hilbert(m_examp_bp(:, k))).^2, 2), 'b-', 'linewidth', 2)
%     
%     output(:, k) = geomean(abs(hilbert(m_examp_bp(:, k))).^2, 2);
%     axis tight,  box off,  ylim([-1.2, 1.2]),
% end
% 
% subplot(8, 1, 8), plot(sum(output, 2)), axis tight
% 
% %% EXTRACT BROADBAND
% 
% % Here, I want to extract broadband using different ways and compare the
% % results.
% bb = {};
% 
% bands = {[60, 200], 20};
% srate = 1000;
% 
% for method = 1 : 5
%    bb{method} = extractBroadband(raw.ts, srate, method, bands);
% end
% 
% % EPOCH THE BROADBAND EXTRACTED -------------------------------------
% bb_epoch = {};
% 
% for method = 1 : 5 % number of broadband extracting method
%    for k = 1 : 71 % number of electrodes
%       tmp = bb{method}(:, k);
%       bb_epoch{method}(:, :, k) = tmp(raw.onsets);
%    end
% end
% 
% %% VISUALIZE THE OUTCOME OF DIFFERENT METHODS AND COMPARE
% 
% % eqaulize the baseline:
% norm_base = @(x) x - mean(x(1 : 200));
% norm_max  = @(x) x./max(x);
% 
% figure (5), clf, color_order = {'g', 'c', 'b', 'k', 'r'};
% 
% for k = 1 : 71
%     subplot(8, 10, k)
%     for method = 1 : 5
%         to_plot = norm_max(norm_base(mean(bb_epoch{method}(:, :, k))));
%         plot(to_plot, '-', 'color', color_order{method}), hold on
%     end
%     axis tight, box off
% end
