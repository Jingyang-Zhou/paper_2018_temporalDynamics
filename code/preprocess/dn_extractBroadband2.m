% dn_extractBroadband version 2

% WHAT IS STORED IN THE RAW DATA ------------------------------------------

% The ECoG experiment consists of 
% (1) 7 stimulus types, each type with 30 repeats, so a total of 210 trials
% (2) each trials consists of 1204 sampling point (time course, 1000Hz)
% (3) a total of 71 good channels and their correpsonding labels
% (4) an extracted broadband matrix, which we will compute here

%% LOAD, MAKE AND SAVE RAW DATA

fName = 'dn_rawData.mat';
dtLoc = fullfile(dn_ECoG_RootPath, 'data');
a     = load(fullfile(dtLoc, fName));
raw   = a.raw;

% ALTERNATIVELY, WE CA RE-GET RAW DATA ------------------------------------
% raw = [];
% [raw.ts, raw.stimNames, raw.imNames, raw.goodLabels, raw.goodChannels,raw.idx] = dn_getRawData;

%% PRE-DEFINED VARIABLES AND FUNCTIONS

% USEFUL FUNCTIONS --------------------------------------------------------
whiten   = @(x) (x - mean(x(:)))./ diff(prctile(x, [.25 .75]));
normMax  = @(x) x./max(x);
normBase = @(x) x-mean(x(1 : 200));

% PREDEFINED VARIABLES ----------------------------------------------------
srate    = 1000;
nElec    = size(raw.ts, 2); % number of electrodes

% HILBERT TRANSFORMS ------------------------------------------------------
stand_hilbert = @(x) abs(hilbert(x));

%% TIME-FREQUENCY ANALYSIS

% In order to decide which bands to use for extracting the broadand signal,
% we are going to do a spectrogram analysis here to visualize how frequency
% componenets changes over time for each electrode.

% 

%% BAND-PASS FILTER THE TIME COURSE IN EACH ELECTRODE

% From Hermes 2015 Figure 1, the example electrode shows th