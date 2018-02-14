% dn_extractBroadband.m
%
% STEPS -------------------------------------------------------------------
% (1) Compute the spectrogram, and decided using two ranges of frequency
%     bands for broadband extraction.
% (2) Band-pass filter using two band ranges ([70, 210], [110, 210]),
%     avoiding the line noise at 16, 120, and 180 Hz.
% (3) Compared two ways of broadband computation, and settled on the
%     simpler way.
%
% THE RAW DATA ------------------------------------------------------------
% The ECoG experiment consists of
% (1) 7 stimulus types, each type with 30 repeats, so a total of 210 trials
% (2) each trials consists of 1204 sampling point (time course, 1000Hz)
% (3) a total of 71 good channels and their correpsonding labels
% (4) an extracted broadband matrix, which we will compute here

%% SAVE KNOBs

saveFigure = 1;
saveData   = 1;

%% LOAD, MAKE AND SAVE RAW DATA

fName = 'dn_rawData.mat';
dtLoc = fullfile(dn_ECoG_RootPath, 'data');
a     = load(fullfile(dtLoc, fName));
raw   = a.raw;

% ALTERNATIVELY, WE CAN RE-GET THE RAW DATA --------------------------------
% raw = [];
% [raw.ts, raw.stimNames, raw.imNames, raw.goodLabels, raw.goodChannels,raw.idx] = dn_getRawData;

%% PRE-DEFINED VARIABLES AND FUNCTIONS

% USEFUL FUNCTIONS --------------------------------------------------------
whiten   = @(x) (x - mean(x(:)))./ diff(prctile(x, [.25 .75]));
normMax  = @(x) x./max(x);
normBase = @(x) x - mean(x(1 : 200));
normBM   = @(x) normMax(normBase(x));

% PREDEFINED VARIABLES ----------------------------------------------------
srate = 1000;
nElec = size(raw.ts, 2); % number of electrodes

% DERIVED VARIABLES -------------------------------------------------------
ts  = raw.ts;
idx = raw.idx;
T   = size(idx, 2)/srate;
t   = 1/srate : 1/srate : T;

% HILBERT TRANSFORMS ------------------------------------------------------
stand_hilbert  = @(x) abs(hilbert(x));
whiten_hilbert = @(x) abs(hilbert(whiten(x)));

%% TIME-FREQUENCY ANALYSIS - SPECTROGRAM

% In order to decide which bands to use for extracting the broadand signal,
% we are going to do a spectrogram analysis here to visualize how frequency
% componenets changes over time for each electrode.

dn_spectrogram_spectra('all', 100); % [may take a little while]

%% BAND-PASS FILTER THE TIME COURSE IN EACH ELECTRODE
bands = {}; bp = {};

% Here we are going to use two band range to extract broadband, [70, 210]
% and [110, 210].

% Wider range
bands{1} = [70 90; 90 110; 110 130; 130 150; 150 170; 190 210];
% narrower range
bands{2} = [110 130; 130 150; 150 170; 190 210];

% BANDPASS FILTER THE ORIGINAL DATA ---------------------------------------
for iband = 1 : length(bands), bp{iband} = bandpassFilter(ts, srate, bands{iband}); end

%% EXTRACT BROADBAND USING THE BANDPASS FILTERED DATA

bb_1 = {}; bb_2 = {}; bb_rs1 = {}; bb_rs2 = {};

% There are multiple ways to extract broadband, here we compare between two ways:
% (1) geomean(asb(hilbert(bp)).^2)
% (2) geomean(abs(hilbert(whiten(bp))).^2)

for iElec = 1 : nElec
    for iband = 1 : length(bands)
        bp_elec = squeeze(bp{iband}(:, iElec, :));
        
        bb_elec1 = geomean(stand_hilbert(bp_elec).^2, 2);
        bb_elec2 = geomean(whiten_hilbert(bp_elec).^2, 2);
        
        % RESHAPE THE EXTRACTED BROADBAND DATA ----------------------------
        bb_rs1{iband}(:, :, iElec) = bb_elec1(idx);
        bb_rs2{iband}(:, :, iElec) = bb_elec2(idx);
    end
end
%% VISUALIZE EXTRACTED BROADBAND AND DECIDE WHICH WAY (OF BROADBAND EXTRACTION) TO USE

% AVERAGE OVER TRIALS BEFORE PLOTTING
for iband = 1 : 2
    mbb_rs1{iband} = squeeze(mean(bb_rs1{iband}));
    mbb_rs2{iband} = squeeze(mean(bb_rs2{iband}));
end

% This figure compares two ways of computing broadband, and each way with
% two band widths. The conclusion is that the difference between the two
% ways of computing broadband is not that big, so we will just pick the
% simpler way (the first way). For other analyses (model fitting and
% parameter estimating), we will do it for both frequency ranges ([70-200],
% [100, 200]).

fg = figure (1); clf
for iElec = 1 : nElec
    subplot(8, 10, iElec)
    plot(t, normBM(mbb_rs1{1}(:, iElec)), 'r-'), hold on
    plot(t, normBM(mbb_rs1{2}(:, iElec)), 'b-'),
    plot(t, 1 + normBM(mbb_rs2{1}(:, iElec)), 'm-'),
    plot(t, 1 + normBM(mbb_rs2{2}(:, iElec)), 'k-'),
    xlim([0, T]), ylim([-0.5, 2]), box off
end
fg.Position = [1, 2000, 2000, 2000];

% COMPARE SPECTROGRAM AND EXTRACTED BROADBAND TIME SERIES -----------------
figure (100),
for iElec = 1 : nElec
    subplot(8, 10, iElec)
    plot(t(201 : 600)-0.13, normBM(mbb_rs2{1}(201 : 600, iElec)).*50 + 50, 'k-', 'linewidth', 3),
    plot(t(201 : 600)-0.13, normBM(mbb_rs2{2}(201 : 600, iElec)).*50 + 50, 'w-', 'linewidth', 3)
end

%% SAVE FIGURE

fNm1    = 'pre_extractBB_methods_bandRng';
fNm2    = 'pre_spectrogram_ts_allstim';
saveLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');

if saveFigure,
    printnice(1, 0, saveLoc, fNm1);
    printnice(100, 0, saveLoc, fNm2);
end

%% SAVE DATA

saveLoc = fullfile(dn_ECoG_RootPath, 'data');
fName   = 'dn_preprocessedData.mat';

if saveData
    dt.ecog     = [];
    dt.ecog.labels = raw.goodLabels;
    dt.ecog.chans  = raw.goodChannels;
    dt.ecog.bb     = bb_rs1;
    dt.ecog.stimNm = raw.stimNames;
    dt.ecog.imNm   = raw.imNames;
    save(fullfile(saveLoc, fName), 'dt')
end
