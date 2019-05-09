% dn_analyzeAdditionalSubject

%% USEFUL FUNCTIONS

whiten         = @(x) (x - mean(x(:)))./ diff(prctile(x, [.25 .75]));
whiten_hilbert = @(x) abs(hilbert(whiten(x))); 

%% LOAD STIMULUS

subjDir   = fullfile('/Volumes', 'server', 'Projects', 'ECoG', 'ecog24-cl', 'ecog', 'Retinotopy');
stimPth   = fullfile(subjDir, 'Stimulus');

stimFiles = dir(fullfile(stimPth, 'SOC_*.mat'));
a         = load(fullfile(stimPth, stimFiles(1).name));

%% ADD ECOG PATH

vista_ecogPth = '~/matlab/git/Vistaproj_ECoG_pRF/';
addpath(genpath(vista_ecogPth))

info     = ecogCreateSession('cl');
timing   = info.timing;
t_start  = timing.starts.soc;
t_end    = timing.ends.soc;

srate    = info.Fs.tseries;

goodElec = info.goodElectrodes;
goodElec = goodElec(goodElec<=118);

%% LOAD AND COMBINE DATA

% SOC scans correspond to file number 15, 16, 17.

data    = {};
dataPth = fullfile(subjDir, 'CARData');
runNums = [15, 16, 17];

for iRun = 1 : 3
    folderPth = fullfile(dataPth, sprintf('S15_86_CL_%d', runNums(iRun)));
    for k = 1 : length(goodElec)
        dtNm = dir(fullfile(folderPth, sprintf('CAR*_%d.mat', goodElec(k))));
        b    = load(fullfile(folderPth, dtNm.name));
        data{iRun}(k, :) = b.wave;
    end
end

%% VISUALIZE IMAGES

% SELECT FULLFIELD/FULLFIELD NOISE STIMULI --------------------------------
fullfield_idx      = [10, 29, 39 : 86];
fullfield_noiseIdx = [10, 29, 55 : 68, 74 : 78, 83 : 86];

figure (1), clf, colormap gray
for k = 1 : 87
    subplot(9, 10, k)
    imagesc(squeeze(a.stimulus.images(:, :, k))), axis off
    if ismember(k, fullfield_noiseIdx), title(k); end
end

%% EXTRACT BROADBAND

raw = [];

%% BAND-PASS FILTER HTE TIME COURSE FOR EACH RUN

bands = [70 80; 80 90; 90 100; 100 110; 130 140; 140 150; 150 160; 160 170; 190 200; 200 210];

% BANDPASS FILTER THE ORIGINAL DATA ---------------------------------------
for irun = 1 : 3
    bp{irun} = bandpassFilter(data{irun}', srate, bands); % time course x n electrodes x 10 bands
end

raw.data = data;
raw.elec = goodElec;
raw.bp   = bp;

%% EXTRACT BROADBAND USING THE BANDPASS FILTERED DATA

bb_elec = {};

for irun = 1 : 3
    for iElec = 1 : length(goodElec)
       bp_elec = squeeze(bp{irun}(:, iElec, :)); 
       bb_elec{irun}(iElec, :) = sum(whiten_hilbert(bp_elec).^2, 2);
    end
end
