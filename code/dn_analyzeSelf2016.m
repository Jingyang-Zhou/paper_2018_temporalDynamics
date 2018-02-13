% Analyze Self multi-unit data

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

%% Dependencies

% Jon's broadband tutorial (or any specific function that I am using)

%% useful functions

normMax_range = @(x, range) (x - mean(x(range)))./max(x - mean(x(range)));
normMax       = @(x) x./max(x);
synapseFunc   = @(x) zeromean(2 * rand(x,1)-1);

%% (1) Load data

bbPth   = '/Users/winawerlab/Google Drive/broadband tutorial';
dataLoc = '/Volumes/server/Projects/Temporal_integration/data/ECoG/';
dataNm  = 'CONTEXT.mat';

a = load(fullfile(dataLoc, dataNm));
addpath(genpath(bbPth))

% a.DETS: [electrode, session number, condition]
% data: 930 samples in length, 400-900 stimulus on

range = 2 : size(a.MUA, 2); % this is the range of data that we will be analyzing

% get electrode and session index
elecIdx = a.DETS(:, 1);
sessIdx = a.DETS(:, 2);

mua{1} = a.MUA(elecIdx == 6, range);
mua{2} = a.MUA(elecIdx == 7, range);

lfp{1} = a.LFP(elecIdx == 6, range);
lfp{2} = a.LFP(elecIdx == 7, range);

m_mua = [];
m_lfp = [];

for k = 1 : 2
    m_mua(k, :) = normMax_range(mean(mua{k}), 1 : 400);
    m_lfp(k, :) = normMax_range(mean(lfp{k}), 1 : 400);
end

%% PRE-DEFINED PARAMETERS

% experimental parameters -------------------------------------------------
lineNoise = 50; % line noise is at 50 Hz
n_elec    = 2;
srate     = a.FsD;
dt        = 1/srate;
t         = dt : dt : size(mua{1}, 2);
t         = t(range);
f_lfpRng  = [1, 200];


n_range = [100 : length(m_mua)-50];

% FOR PART 1: process LFP signal ------------------------------------------
epochLth = 0.4; % second
iBase    = t <= epochLth + dt;
stimOn   = 0.4;
iStim    = t > (stimOn + 0.05) & t - (stimOn + 0.05) < epochLth;
nf       = sum(iStim);
f        = [0 : nf - 1] ./ epochLth;
t_plot   = t - stimOn;

% FOR PART 1.2 : for computing the spectrogram ----------------------------
spectro_prms = [];
spectro_prms.pad      = -1;
spectro_prms.tapers   = [3 5];
spectro_prms.fpass    = [0 200];
spectro_prms.Fs       = a.FsD;
spectro_prms.trialave = 1;
movingwin             = [.100 .01]; % length of the moving window and step size


% FOR PART 1.3 : for computing LFP broadband ------------------------------
bands     = {[100, 170], 20}; % {[lower bound,  upper bound], window sz}

% for individual bands:
bands_indi = {}; idx = 0;

for k = f_lfpRng(1) : 30 : f_lfpRng(2)
    idx = idx + 1; bands_indi{idx} = {[k, k + 30], 10};
    bands_indi_str{idx} = sprintf('[%d, %d]', k, k + 30);
end

% FOR PART 2 : dendritic model fit ----------------------------------------

% for model fit to MUA data -----------------------------------------------
dn_init  = [0.04, 0.1, 1.8, 0.05, 0.01, 1]; % 'tau1',  'tau2', 'n', 'sigma', 'shift', 'scale'
dn_lb    = [0,0,0,0,0,0];
dn_ub    = [1, 1, 10, 1, 1, 10];

% create stimulus
stim = zeros(1, length(lfp{1}));
stim(400 : 900) = 1;
stim = stim(n_range);

% for the dendritic current analysis --------------------------------------
spikeRateScale = 0.1;
nSynapse       = 2000;
nTrials        = 100;
restSpike      = 0.1; % final resting spike is the product of restSpike and spikeRateScale

% two time constants:

alpha   = 0.1; % % time constant for dendritic integration, from Miller et al. (0.1)
tau     = 0.0023; % time constant for post-synaptic current, from Miller et al. (0.0023)

nt      = length(n_range);
nTrials = 50;

% for dendritic current fit -----------------------------------------------
nRounds = 100;
nTrials = 80;
init    = [tau, alpha];
lb      = [10^(-5), 10^(-5)];
ub      = [1, 1];

%% PART 1: SIGNAL PROCESSING

%% PART 1.1 : process LFP signal

f_lfpBase = {}; f_lfpStim = {}; notched_lfp = {}; mf_lfpBase = []; mf_lfpStim = [];

% fft the baseline portion and the stimulus-on portion of the time series
for k = 1 : n_elec
    for k1 = 1 : size(lfp{k}, 1)
        % apply notch filter to LFP
        notched_lfp{k}(k1, :) = dn_ecogNotch(lfp{k}(k1, :)', srate, lineNoise);
        
        % compute frequency in the baseline and the stimulus on range
        f_lfpBase{k}(k1, :) = abs(fft(notched_lfp{k}(k1, iBase)));
        f_lfpStim{k}(k1, :) = abs(fft(notched_lfp{k}(k1, iStim)));
    end
    % average frequency
    mf_lfpBase(k, :) = mean(f_lfpBase{k});
    mf_lfpStim(k, :) = mean(f_lfpStim{k});
end

disp('PART 1.1 : process LFP signal')

%% PART 1.2 : Compute spectrogram (from 100 Hz and up)

spectro = {}; spectrobase = {}; m_spectro = []; m_spectrobase = [];

for k = 1 : n_elec
    for k1 = 1 : size(lfp{k}, 1)
        [spectrobase{k}(k1, :,:), t_spectro, f_spectro] = mtspecgramc(notched_lfp{k}(k1, iBase), movingwin, spectro_prms);
        [spectro{k}(k1, :,:), t_spectro, f_spectro] = mtspecgramc(notched_lfp{k}(k1, iStim), movingwin, spectro_prms);
    end
    
    % take the trial mean and normalize spectrogram to the base
    m_spectrobase(k, :, :) = mean(spectrobase{k});
    m_spectro(k, :, :)     = mean(spectro{k})./m_spectrobase(k, :, :);
end

disp('PART 1.2 : Compute spectrogram')

%% VISUALIZATION: spectrum and spectrogram

figure (2), clf, subplot(3, 3, 1), linetype = {'-', '-.'};

t_spectro = t_plot(find(iStim, 1)) + t_spectro - t_spectro(1);

for k = 1 : 2
    loglog(f, mf_lfpBase(k, :), sprintf('k%s', linetype{k}), 'linewidth', 3), hold on
    loglog(f, mf_lfpStim(k, :), sprintf('r%s', linetype{k}), 'linewidth', 3),
end
patch([20, 186, 186, 20], [60, 60, 10^3.8, 10^3.8], 'y', 'facealpha', 0.2, 'edgecolor', 'r'),
patch([37, 60, 60, 37], [60, 60, 10^3.8, 10^3.8], 'b', 'facealpha', 0.15, 'edgecolor', 'k'),
% use box to indicate the broaband range and the gamma oscillation range
legend('baseline, elec 6', 'stim on, elec 6', 'baseline, elec 7', 'stim on, elec 7', ...
    'location', 'southwest')
xlabel('frequency'), ylabel('amplitude'), title('LFP spectrum')
axis tight, xlim([0, nf/2]), box off, set(gca, 'fontsize', 13)

%  PART 1.2: visualize spectrogram
for k = 1 : 2
    subplot(3, 3, k + 1), imagesc(t_spectro, f_spectro, log10(squeeze(m_spectro(k, :, :)))', [-.7, .7]), hold on
    patch([t_spectro(1), t_spectro(end), t_spectro(end), t_spectro(1)], [100, 100, max(f), max(f)], 'r', 'facealpha', 0, 'edgecolor', 'r', 'linewidth', 3),
    box off, axis xy, set(gca, 'fontSize', 14),
    xlabel('time (s)'), ylabel('frequency (Hz)'), title(sprintf('electrode %d (stim on)', k + 5))
end

%% PART 1.3.1: compute LFP broadband

bb  = {}; mbb = []; 

% option "4": geomean(abs(hilbert(whiten(bp(x)))).^2)'
for k = 1 : n_elec
    bb{k}(:, :)     = extractBroadband(notched_lfp{k}', srate, 4, bands);
    mbb(k, :)       = normMax_range(median(bb{k}(n_range, :)'), 1 : 300);
end

%% PART 1.3.2: compute LFP broadband (every 30 Hz)

bb_indi = {}; mbb_indi = [];

for k = 1 : n_elec
    for k1 = 1 : length(bands_indi)
        bb_indi{k}{k1} = extractBroadband(notched_lfp{k}', srate, 4, bands_indi{k1});
        mbb_indi(k, k1, :) = normMax_range(median(bb_indi{k}{k1}(n_range, :)'), 1 : 300);
    end
end

%% PART 1.3.3: fit normalization model to different LFP bands

model = [];

for k = 1 : n_elec
    for k1 = 1 : length(bands_indi)
        model.LFP_indiPrm(k, k1, :) = fminsearchbnd(@(x) dn_computeFineFit(x, squeeze(mbb_indi(k, k1, :))',...
            stim, t(n_range), 'uniphasic'), dn_init, dn_lb, dn_ub);
    end
end

%% VISUALIZE : BROADBAND

for k = 1 : n_elec
    figure (2), subplot(3, 3, k + 4)
    plot(t_plot(n_range), m_mua(k, n_range), 'k-'), hold on
    plot(t_plot(n_range), mbb(k, :),'b-', 'linewidth', 4),

    set(gca, 'ytick', [0, 0.5, 1])
    axis tight, box off, xlabel('time (s)'), ylabel('normalized amplitude'), title(sprintf('Electrode %d time series', k + 5))
    set(gca, 'fontsize', 12), legend('MUA', 'LFP broadband', 'LFP broadband (low)', 'Location', 'southeast')
end

figure (2),
for k = 1 : n_elec
    
    subplot(3, 3, k + 7); cla,  set(gca, 'ColorOrder', copper(length(bands_indi))),hold on
    stem(t_plot(n_range), m_mua(k, n_range), 'r-', 'markersize', 1), hold on
    plot(t_plot(n_range), squeeze(mbb_indi(k, 2 : end , :))', 'linewidth', 2), axis tight
    
    legend(bands_indi_str{2 : end}, 'location', 'southwest'), set(gca, 'fontsize', 12)
end

%%

%% PART 2 USE MUA TO PREDICT LFP

%% PART 2.0 : generate model predictions to MUA time course

model = [];

for k = 1 : n_elec
    model.MUAprm(k, :) = fminsearchbnd(@(x) dn_computeFineFit(x, m_mua(k, n_range), stim(n_range), t(n_range), 'uniphasic'), dn_init, dn_lb, dn_ub);
    [~,~, model.MUAprd(k, :)] = dn_computeFineFit(model.MUAprm(k, :), m_mua(k, n_range), stim(n_range), t(n_range), 'uniphasic');
end

%% PART 2.1 : visualize MUA model fit

figure (1), clf
for k = 1 : n_elec
    subplot(1, 2, k)
    plot(t(n_range), m_mua(k, n_range), 'k-'), hold on
    plot(t(n_range), model.MUAprd(k, :), 'r-', 'linewidth', 2)
    xlabel('time (s)'), ylabel('normalized amplitude'), set(gca, 'ytick', [0, 0.5, 1])
    title(sprintf('electrode %d', k + 5)), set(gca, 'fontsize', 15),  box off, axis tight, axis square
end

%% PART 2.1 : variable parameters


%% PART 2.1 simplified way to generate dendritic current

% (1) generate a random variable that takes values only from [-1, 1],
% representing the spikes from excitatory and inhibitory neurons
mua_type     = [];
I            = {};
pred_bb_simp = {};

for k = 1 : n_elec
    for k1 = 1 : nTrials
        neuronType = rand(1, nt);
        neuronType(neuronType > 0.5)  = 1;
        neuronType(neuronType <= 0.5) = -1;
        mua_type{k}(k1, :) = m_mua(k, n_range).*neuronType;
    end
end

% (2) compute dendritic current
for k = 1 : n_elec
    I{k} = dn_computeDendriticCurrent(mua_type{k}, tau, alpha, dt);
end

% (3) extract broadband
for k = 1 : n_elec
    for k1 = 1 : nTrials
        pred_bb_simp{k}(k1, :) = extractBroadband(I{k}(k1, :)', srate, 4, bands);
    end
end

figure (5), clf
subplot(1, 2, 1), stem(t(n_range), normMax_range(m_mua(1, n_range), [1 : 300]), 'k-', 'markersize', 1), hold on
plot(t(n_range), normMax_range(median(pred_bb_simp{1}), [1 : 300]), 'r-', 'linewidth', 3),
plot(t(n_range), mbb(1, n_range), 'b-', 'linewidth', 3), box off

subplot(1, 2, 2), stem(t(n_range), normMax_range(m_mua(2, n_range), [1 : 300]), 'k-', 'markersize', 1), hold on
plot(t(n_range), normMax_range(median(pred_bb_simp{2}), [1 : 300]), 'r-', 'linewidth', 3),
plot(t(n_range), mbb(2, n_range), 'b-', 'linewidth', 3), box off


%% PART 2.2 more complex way to generate dendritic current
param = [];

% create Poisson spike rate
spikeRate   = {};
totalSpikes = {};

for k = 1 : n_elec
    model.MUAprdRest(k, :) = normMax(model.MUAprd(k, :) + restSpike).*spikeRateScale;
end

% create Poisson rate
for k = 1 : n_elec
    spikeRate{k} = repmat(model.MUAprdRest(k, :)', 1, nSynapse);
end

for iRound = 1 : nRounds
    disp(sprintf('round %d', iRound))
    for k = 1 : n_elec
        % compute total Spikes
        totalSpikes{iRound}(k, :, :) = dn_generatePoissonSpikes(spikeRate{k}, nTrials, nSynapse);
        
        % generate dendritic current
        [param(iRound, k, :), ~, exitflag(iRound, k)] = fminsearchbnd(@(x)dn_fit_dendriticCurrent(squeeze(totalSpikes{iRound}(k, :, :))', mbb(k, n_range), x, dt, srate), init, lb, ub);
    end
end

%% PART 2.3 : post-processing

% alpha is the constant for dendritic integration, and tau is the
% post-synaptice time constant
I        = {};
fit_bb   = {};
m_fit_bb = [];

% compute median parameters
for k = 1 : n_elec
    m_alpha(k) = median(param(:, k, 1));
    s_alpha(k, :) = prctile(param(:, k, 1), [25, 75]);
    
    m_tau(k) = median(param(:, k, 2));
    s_tau(k, :) = prctile(param(:, k, 2), [25, 75]);
end

% compute predicted LFP broadband
for iRound = 1 : nRounds
    for k = 1 : n_elec
        I{iRound}(k, :, :) = dn_computeDendriticCurrent(squeeze(totalSpikes{iRound}(k, :, :)), param(iRound, k, 1),  param(iRound, k, 2), dt);
        for k1 = 1 : nTrials
            fit_bb{iRound}(k, k1, :) = extractBroadband(squeeze(I{iRound}(k, k1, :))', srate, 4, bands);
        end
    end
end

% compute the mean fitted broadband
for iRound = 1: nRounds
    m_fit_bb(iRound, 1, :) = median(fit_bb{iRound}(1, :, :), 2);
    m_fit_bb(iRound, 2, :) = median(fit_bb{iRound}(2, :, :), 2);
end

%% PART 2.3 : visualize parameter distributions

figure (3), clf

% first plot distribution of the two parameters
subplot(2, 3, 1), location = [1, 1.2];
for k = 1 : 2
    yyaxis left
    plot(location(k), m_alpha(k), 'bo', 'markerfacecolor', 'b', 'markersize', 8), hold on
    plot([location(k), location(k)], [s_alpha(k, 1), s_alpha(k, 2)], 'b-', 'linewidth', 3), ylim([2 * 10^(-3), 3 * 10^(-3)])
    
    yyaxis right
    plot(location(k) + 0.5, m_tau(k), 'ro', 'markerfacecolor', 'r', 'markersize', 8), hold on
    plot([location(k) + 0.5, location(k) + 0.5], [s_tau(k, 1), s_tau(k, 2)], 'r-', 'linewidth', 3)
end
xlim([0.9, 1.8]), box off, set(gca, 'fontsize', 14, 'xtick', [1, 1.2, 1.5, 1.7], 'xticklabel', {'alpha', 'alpha', 'tau', 'tau'})
title('Time constants')

% second plot predicted LFP versus derived LFP
for k = 1 : 2
    subplot(2, 3, k + 1)
    %stem(t(n_range), m_mua(k, n_range), 'k-', 'markersize', 1), hold on
    plot(t(n_range(1 : cutoff)), model.MUAprd(k, (1 : cutoff)), 'k-'), hold on
    plot(t(n_range(1 : cutoff)), normMax_range(median(bb{k}(:, n_range(1 : cutoff))), [1 : 300]), 'b-', 'linewidth', 3),
    plot(t(n_range(1 : cutoff)), normMax_range(squeeze(median(m_fit_bb(:, k, (1 : cutoff)), 1))', [1 : 300]), 'r-', 'linewidth', 3);
end

%%

toplot = normMax_range(mean(fit_bb{1}(k, :, :), 2), [1 : 300]);

figure (5), clf
for k = 1 : 2
    subplot(1, 2, k), 
    stem(t(n_range(1 : cutoff)), m_mua(k, n_range(1 : cutoff)), 'k-', 'markersize', 1), hold on
    plot(t(n_range(1 : cutoff)),  model.MUAprd(k, (1 : cutoff)), 'linewidth', 4)
end
