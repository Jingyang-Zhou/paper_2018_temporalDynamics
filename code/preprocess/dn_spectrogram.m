% dn_spectrogram

% INPUTS ------------------------------------------------------------------

% ts  : n time points x number of electrodes 
% idx : re-shape the time series into n trials x time course per trial

% OUTPUTS -----------------------------------------------------------------

%% EXAMPLE

fName = 'dn_rawData.mat';
dtLoc = fullfile(dn_ECoG_RootPath, 'data');
a     = load(fullfile(dtLoc, fName));
raw   = a.raw;

ts  = raw.ts;
idx = raw.idx;

%% DERIVED PARAMETERS 

n_elec = size(ts, 2); % number of electrodes

%% RE-SHAPE THE DATA INTO N_TRIALS X X_TIMEPOINTS_PER_TRIAL
ts_rs = [];

% re-shape the time course
for k = 1 : n_elec
    ts_elec = ts(:, k);
    ts_rs(:, :, k) = ts_elec(idx);
end

%% DEFINE TIME-FREQUENCY MULTITAPER SETTINGS

% same setting as in Hermes et al. 2015, except that we average data across
% 3 stimulus types

movingwin=[.200 .05];
params.pad=-1;
params.tapers=[3 5];
params.fpass=[0 200];
params.Fs=srate;
params.trialave=1;

%% CAUCULATE THE SPECTROGRAM


