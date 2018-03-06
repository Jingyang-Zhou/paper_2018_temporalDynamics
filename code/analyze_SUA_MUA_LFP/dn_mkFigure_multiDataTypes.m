
% TO DO: WHEN BOOTSTRAPPING, TAKE THE NUMBER OF CELLS IN EACH TYPE INTO
% ACCOUNT

%%

% dn_Alternative measurement

% In this file, we fit the DN model to three data types:
% (1) average of single unit data (Albrect and Geisler 2002?)
% (2) MUA data, averaged over 2 electrodes (Self et al. 2017)
% (3) LFP data, averaged over 2 electrodes (Self et al. 2017)

%% MULTI-UNIT DATA

% en_extractBroadband multi-Unit
% DATA DESCRIPTION -------------------------------------------------

% two electrodes in area V2/V3

% data format for MUA and LFP: 2D matrix: [trials, samples]

% The raw signals from E6 and E7 (the two electrodes) were re-referenced to
% the average of the (nonspiking) electrodes E0-E4 (E5 was excluded from
% the average due to high noise levels). From the re-referenced signal they
% created three signals: the local field potential (LFP), the envelop of
% multi-unit activity (MUA) and the thresholded multi-unit activity (MUAt).
% The LFP was created by first down-sampling to 930 Hz, then band-pass
% filtering the resulting signal between 1 Hz and 200 Hz using a second
% order, zero-phase Butterworth filter.

% Line-noise was removed by fitting a 50 Hz sine-wave to each individual
% trial, then subtracting it.

% They measured MUA by band-passing the raw signal between 500 Hz and 5 kHz
% to isolate high-frequency (spiking) activity. This filtered signal was
% rectified, down-sampled to 930 Hz and low-pass filtered (<200 Hz) to
% measure the envelope of the spiking activity. They also generated MUAt by
% thresholding the band-passed signal.

% (Self et al. 2016, page 17/26)

% stimulus: They measured RFs by flashing a small black-and-white
% checkerboard pattern(1 x 1 deg. in size), check size (0.33 deg.) at every
% point of an 11 x 11 deg. grid.

% contextual modulation experiment: stationary gratings. The phase of the
% central gratin was randomly chosen on each trial from a uniform
% distribution ranging from -pi to pi. They used two different, orthogonal
% orientations for the central grating, 60 and 150 deg. They presented an
% equal number of trials with the two orientations so that the average
% stimulation of the RF was identical for all conditions. The stimulus
% durations was 500 ms with an inter-trial interval of 500 ms.

% Contents of the DETS variables:
% Receptive field mapping: [electrode, check position index, X-position, Y-position]
% positions are given in degrees

% Contextual effects: [electrode, session number, condition]

% link: https://osf.io/2euy6/

% (not applied). RF tuning properties wewre measured using drifting
% sine-wave gratings placed at a location that activated neurons at both E6
% and E7 (10.3 DEG. ECCENTRICITY, -14 DEG. ANGLE FROM HORIZONTAL MERIDIAN)

%% USEFUL FUNCTIONS

% HILBERT TRANSFORMS ------------------------------------------------------
stand_hilbert  = @(x) abs(hilbert(x));

%% LOAD DATA AND MODEL PARAMS(CONTEXTUAL EFFECT)

% LOAD SELF DATA ----------------------------------------------------------
dtNm  = 'Self_CONTEXT.mat';
dtLoc = fullfile(dn_ECoG_RootPath, 'data');
a     = load(fullfile(dtLoc, dtNm));

% LOAD MODEL PARAMS -------------------------------------------------------
prm_fNm  = 'dn_params.mat';
b        = load(fullfile(dtLoc, prm_fNm));

% LOAD SINGLE UNIT DATA ---------------------------------------------------
s_dtNm = 'figure4Data.xlsx';
c      = xlsread(fullfile(dtLoc, s_dtNm));

% LOAD ECOG BROADBAND DATA ------------------------------------------------
e_dtNm = 'dn_params.mat';
d      = load(fullfile(dtLoc, e_dtNm));

%% PRE-DEFINED OR DERIVED VARIABLES

mua = []; lfp = [];

% ECOG PARAMETERS ---------------------------------------------------------
ecog = [];
ecog.dt  = normMax(d.prm.ecog.dn.bs_bbts(:, 1 : 100));
ecog.prd = normMax(d.prm.ecog.dn.prd (1 : 100, :)');
ecog.t   = [1 : size(ecog.dt)]./1000;
ecog.stim = zeros(1, length(ecog.t));
ecog.stim(201 : 700) = 1;

% EXTRACTED PARAMETERS ----------------------------------------------------
mua.dt  = a.MUA;  % MUA data, data format: trials x samples in time
lfp.dt  = a.LFP;  % LFP data, data format: trials x samples in time
dtInfo  = a.DETS; % [electrode, session number, condition]
srate   = a.FsD;
nTrials = size(mua.dt, 1);

% DERIVED PARAMETERS ------------------------------------------------------
dt = 1/srate;

% make time points:
T = size(mua.dt, 2) .* dt;
t = dt : dt : T;

% make stimulus:
stim = zeros(1, length(t));
stimStart = 0.4; stimEnd = 0.9;
stim(t > stimStart & t <= stimEnd) = 1; % 500 ms stimulus


trimIdx = 251 : 930; % for visualization, get rid of some of the baseline repsonse period
t       = t(trimIdx) - 0.25;
stim    = stim(trimIdx);
mua.dt  = mua.dt(:, trimIdx);

% FOR THE SINGLE UNIT DATA ------------------------------------------------
sua = [];
sua.t  = c(:, 1)./1000;
sua.dt = c(:, 2 : end);
sua.stim = ones(1, length(sua.t));
nCells = size(sua.dt, 2);
nBoots = 100;

% DN MODEL FITTING --------------------------------------------------------
irfType = 'uniphasic';
dn      = b.prm.mua.dn;

%% ANALYZE SUA:

for iBoot = 1 : nBoots
    idx = randi(nCells, [1, nBoots]);
    sua.bsdt(iBoot, :) = normMax(mean(sua.dt(:, idx), 2));
end
sua.mbsdt = mean(sua.bsdt);
sua.sbsdt = std(sua.bsdt);

figure (1), clf
subplot(1, 4, 1), title('SUA (macaque V1)'), hold on
patch([0, 0.2, 0.2, 0], [0, 0, 1, 1], 0.8 * ones(1, 3)),
shadedErrorBar(sua.t, sua.mbsdt, sua.sbsdt, 'b'), xlim([-0.15, 0.3]), ylim([-0.1, 1]), set(gca, 'fontsize', 14)
set(gca, 'xaxislocation', 'origin'), set(gca, 'xtick', [0, 0.2], 'ytick', [0, 1]), axis square
ylabel('normalized response '),

%% BOOTSTRAP THE SUA DATA

% COARSE FIT --------------------------------------------------------------
[sua.seed, sua.seedR2] = dn_gridFit(sua.mbsdt, dn, sua.stim, sua.t, irfType);

% FINE FIT ----------------------------------------------------------------
sua.seed = [sua.seed, 0];
[sua.prm, ~, sua.r2, sua.exitFlg] = dn_fineFit(sua.mbsdt, sua.stim, sua.t, dn, sua.seed, irfType);

% MAKE PREDICTION ---------------------------------------------------------
sua.prm = [sua.prm(1), 0, sua.prm(2 : end), 1];
sua.prd = normMax(dn_DNmodel(sua.prm, sua.stim, sua.t));

% PLOT MODEL FIT ----------------------------------------------------------
figure (1), subplot(1, 4, 1), plot(sua.t, sua.prd, 'r', 'linewidth', 2)


%% ANALYZE MUA:

%% COMPUTE THE AVERAGED MUA

mua.mdt = mean(mua.dt);

% GET RID OF THE BASELINE -------------------------------------------------
baseline = mean(mua.mdt(t <= stimStart - 0.25));
mua.mdt  = mua.mdt - baseline;

% NORMALIZE THE AVERAGED RESPONSE TIME COURSE TO ITS MAX ------------------
mua.mdt = normMax(mua.mdt')';

%% VISUALIZE THE AVERAGED MUA

figure (1), 
subplot(1, 4, 2), title('MUA (human V2/V3)')
% PLOT THE STIMULUS -------------------------------------------------------
patch([0.15, 0.65, 0.65, 0.15], [0, 0, 1, 1], 0.8 * ones(1, 3)), hold on

plot(t, mua.mdt, 'b-', 'linewidth', 1, 'markersize', 1), axis tight, 
set(gca, 'xaxislocation', 'origin', 'fontsize', 14, 'ytick', [0, 1], 'xtick', [0.15, 0.65], 'xticklabel', [0, 0.5])
axis square, ylim([-0.1, 1])
 
%% FIT DN MODEL TO THE AVERAGED MUA DATA

% COARSE FIT --------------------------------------------------------------
[mua.seed, mua.seedR2] = dn_gridFit(mua.mdt, dn, stim, t, irfType);

% FINE FIT ----------------------------------------------------------------
mua.seed = [mua.seed, 0];
[mua.prm, ~, mua.r2, mua.exitFlg] = dn_fineFit(mua.mdt, stim, t, dn, mua.seed, irfType);

% MAKE PREDICTION ---------------------------------------------------------
mua.prm = [mua.prm(1), 0, mua.prm(2 : end), 1];
mua.prd = normMax(dn_DNmodel(mua.prm, stim, t));

% PLOT MODEL FIT ----------------------------------------------------------
figure (1), subplot(1, 4, 2), plot(t, mua.prd, 'r', 'linewidth', 3)

%% ANALYZE LFP: 

%% STEPS TO EXTRACT BROADBAND

% (1) Because the authors used grating stimulus, need to look at spectra
%     and spectrogram to decide which range of frequency bands contains the
%     gamma oscillation.
% (2) Extract broadband

%% LFP : SPECTROGRAM (CURRENTLY SKIPPING THIS STEP BECAUSE I'VE DONE IT BEFORE) / SPECTRA

fft_dt = []; base_dt = [];
for k = 1 : size(lfp.dt, 1)
    fft_dt(k, :) = abs(fft(lfp.dt(k, 401 : end)));
    base_dt(k, :) = abs(fft(lfp.dt(k, 1 : 400)));
end

nf1 = length(401 : 930);
nf2 = length(1 : 400);

f1 = [0 : nf1 - 1]./nf1 * 1000;
f2 = [0 : nf2 - 1]./nf2 * 1000;

figure (10), clf
loglog(f1, mean(fft_dt), 'k-', 'linewidth', 3), hold on
loglog(f2, mean(base_dt), '-', 'color', 0.7 * ones(1, 3), 'linewidth', 3), xlim([20, 140]), box off
legend('baseline', 'stimulus dependent'), set(gca, 'fontsize', 16)
set(gca, 'xtick', [20, 60, 100, 200])
%% EXTRACT BROADBAND

bands = [75 95; 105 125; 125 145; 155 175]; % line noise is at 50 Hz

% BANDPASS-FILTER THE LFP TIME COURSES ------------------------------------
bp = bandpassFilter(lfp.dt', srate, bands);

% EXTRACT BROADBAND -------------------------------------------------------
for iTrial = 1 : nTrials
    for iband = 1 : size(bands, 1)
        bp_trial = squeeze(bp(:, iTrial, iband));
        bb0(iTrial, iband, :) = geomean(stand_hilbert(bp_trial).^2, 2);
    end
end

bb = squeeze(mean(bb0, 2));

%% PLOT THE LFP BROADBAND DATA

mbb = mean(bb);
mbb = mbb(trimIdx); 
baseline = mean(mbb(t<0.15));
lfp.mbb = normMax(mbb - baseline);

figure (1)
subplot(1, 4, 3), cla, title('LFP (human V2/V3)'), hold on
% PLOT THE STIMULUS -------------------------------------------------------
patch([0.15, 0.65, 0.65, 0.15], [0, 0, 1, 1], 0.8 * ones(1, 3)),
% PLOT THE LFP TIME SERIES ------------------------------------------------
plot(t(1 : end - 36), lfp.mbb(1 : end - 36), 'b-', 'linewidth', 3), axis square
set(gca, 'xaxislocation', 'origin', 'fontsize', 14, 'ytick', [0, 1], 'xtick', [0.15, 0.65], 'xticklabel', [0, 0.5]),
axis tight, ylim([-0.1, 1])

%% FIT THE DN MODEL TO THE LFP BROADBAND

lfp_idx = 1 : length(t) - 35;

% COARSE FIT --------------------------------------------------------------
[lfp.seed, lfp.seedR2] = dn_gridFit(lfp.mbb(lfp_idx), dn, stim(lfp_idx), t(lfp_idx), irfType);

% FINE FIT ----------------------------------------------------------------
lfp.seed = [lfp.seed, 0];
[lfp.prm, ~, lfp.r2, lfp.exitFlg] = dn_fineFit(lfp.mbb(lfp_idx), stim(lfp_idx), t(lfp_idx), dn, lfp.seed, irfType);

% MAKE PREDICTION ---------------------------------------------------------
lfp.prm = [lfp.prm(1), 0, lfp.prm(2 : end), 1];
lfp.prd = normMax(dn_DNmodel(lfp.prm, stim(lfp_idx), t(lfp_idx)));

% PLOT MODEL FIT ----------------------------------------------------------
figure (1), subplot(1, 4, 3), plot(t(lfp_idx), lfp.prd, 'r', 'linewidth', 2)

%% PLOT ECOG DATA AND MODEL FIT

figure (1), subplot(1, 4, 4), cla, title('ECOG broadband (human V1)')
patch([0, 0.5, 0.5, 0], [0, 0, 1, 1], 0.8 * ones(1, 3)), hold on
shadedErrorBar(ecog.t - 0.2, mean(ecog.dt, 2), std(ecog.dt, [], 2), 'b-'), 
plot(ecog.t - 0.2, mean(ecog.prd, 2), 'r-', 'linewidth', 2)
xlim([-0.1, 0.6]), ylim([-0.1, 1]), set(gca, 'xaxislocation', 'origin'), axis square
set(gca, 'fontsize', 14, 'xtick', [0, 0.5], 'ytick', [0, 1])

%% SAVE FIGURE

figLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
figNm  = 'ana_DNFit2modalities';
printnice(1, 0, figLoc, figNm);
