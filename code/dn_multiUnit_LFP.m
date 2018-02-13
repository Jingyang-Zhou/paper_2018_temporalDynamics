% dn_multiUnit_LFP

% TO DO: REMOVE GAMMA OSCILLATION AND COMPUTE THE BROADBAND SIGNAL AGAIN

% QUESTION: IS MUA SERIES IN UNIT OF SPIKES PER 1/SRATE?


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

%% dependencies

bbPth = '/Users/winawerlab/Google Drive/broadband tutorial';

addpath(genpath(bbPth))

%%

dataLoc = '/Volumes/server/Projects/Temporal_integration/data/ECoG/';

normMax = @(x) (x - mean(x(1 : 100)))./max(x - mean(x(1 : 100)));
normMax_noBase = @(x) x./max(x);

%% load Self data (receptive field mapping)


dataNm  = 'RFmapping.mat';
a       = load(fullfile(dataLoc, dataNm));

mua = a.MUA;
lfp = a.LFP;
DETS = a.DETS;

% plot data:
% m_mua = mean(mua);
m_lfp = mean(lfp);

% compute the mean mua responses
m_mua = [];
range = [100: 400];

m_mua(1, :) = mean(mua(DETS(:, 1)==6, range));
m_mua(2, :) = mean(mua(DETS(:, 1) == 7, range));

for k = 1 : size(m_mua, 1)
    m_mua(k, :) =  normMax(m_mua(k, :) - mean(m_mua(k, :)));
end

figure (1), clf
subplot(1, 2, 1)
plot(m_mua(1, :), 'k-'), hold on,
subplot(1, 2, 2)
plot(m_mua(2, :), 'k-'), hold on
axis tight, box off

%% preprocess LFP (receptive field mapping)

% frequency bands for extracting broadband
bands = {[70 170], 20}; % {[lower bound,  upper bound], window sz}

for k = 1 : size(lfp, 1)
    r_trialbb(k, :) = extractBroadband(lfp(k, :), srate, 4, bands);
end

mbb_r = [];

mbb_r(1, :) = normMax(median(r_trialbb(c_DETS(:, 1) == 6, range)));
mbb_r(2, :) = normMax(median(r_trialbb(c_DETS(:, 1) == 7, crange)));

%% 

figure (1),
subplot(1, 2, 1)
plot(mbb_r(1, :), 'b-')
subplot(1, 2, 2)
plot(mbb_r(2, :), 'b-')



%% load Self data (contextual effect)

% need to remove possible Gamma oscillation
% stimulus on : [400, 900] in the original data or [150, 650] on the figure

dataNm1 = 'CONTEXT.mat';
b       = load(fullfile(dataLoc, dataNm1));

c_DETS = b.DETS;
c_mua  = b.MUA;
c_lfp  = b.LFP;
srate  = b.FsD;

dt = 1/srate;

t = dt : dt : size(c_mua, 2)*dt;

cm_mua      = [];
cm_mua_orig = [];

crange = [250 : 930];

stim = zeros(1, length(t));
stim(401 : 900) = 1;

cm_mua(1, :) = normMax(mean(c_mua(c_DETS(:, 1) == 6, crange)));
cm_mua(2, :) = normMax(mean(c_mua(c_DETS(:, 1) == 7, crange)));

cm_mua_orig(1, :) = mean(c_mua(c_DETS(:, 1) == 6, :));
cm_mua_orig(2, :) = mean(c_mua(c_DETS(:, 1) == 7, :));

% for k = 1 : size(cm_mua, 1)
%     cm_mua(k, :) =  normMax(cm_mua(k, :) - mean(cm_mua(k, :)));
% end

figure (2), clf
subplot(2, 2, 1)
plot(t(crange), cm_mua(1, :), 'k-', 'linewidth', 2), hold on, title('electrode 6')

subplot(2, 2, 2)
plot(t(crange), cm_mua(2, :), 'k-', 'linewidth', 2), hold on, title('electrode 7')

%% preprocess LFP

% frequency bands for extracting broadband
bands = {[70 170], 20}; % {[lower bound,  upper bound], window sz}

for k = 1 : size(c_lfp, 1)
    trialbb(k, :) = extractBroadband(c_lfp(k, :), srate, 4, bands);
end

%% plot LFP

mbb = [];

mbb(1, :) = normMax(median(trialbb(c_DETS(:, 1) == 6, crange)));
mbb(2, :) = normMax(median(trialbb(c_DETS(:, 1) == 7, crange)));

figure (2), 
subplot(2, 2, 1), plot(t(crange), mbb(1, :), 'b-', 'linewidth', 3), 
legend('MUA', 'LFP broadband'), xlabel('time (s)'), box off, set(gca, 'fontsize', 12)

subplot(2, 2, 2), plot(t(crange), mbb(2, :), 'b-', 'linewidth', 3), 
legend('MUA', 'LFP broadband'), xlabel('time (s)'), box off, set(gca, 'fontsize', 12)

%% USE MUA TO PREDICT LFP HERE:

%% step 1: Generate model prediction to mean mua time course

init  = [0.07, 0.1, 1.5, 0.02, 0, 1]; % 'tau1',  'tau2', 'n', 'sigma', 'shift', 'scale'
model = [];

for k = 1 : 2
   model.MUAprm(k, :) = fminsearch(@(x) dn_computeFineFit(x, cm_mua(k, :), stim(crange), t(1 : length(crange)), 'uniphasic'), init);
   [~,~, model.MUAprd(k, :)] = dn_computeFineFit(model.MUAprm(k, :), cm_mua(k, :), stim(crange), t(1 : length(crange)), 'uniphasic');
end

figure (2)
for k = 1 : 2
    subplot(2, 2, k), plot(t(crange), model.MUAprd(k, :), 'r--', 'linewidth', 2)
end

%% Use MUA to predict LFP here -- step 1: generate poisson rate

nt        = length(t(crange));
nSynapses = 1000;

restingSpikeRate = 0.1;

for k = 1 : 2
   model.MUAprdRest(k, :) = normMax_noBase(model.MUAprd(k, :) + 0.1);
end

% poisson spike rate is the same as the MUA sequence for now
rate = {};
rate{1} = repmat(model.MUAprdRest(1, :).*0.2, nSynapses, 1)';
rate{2} = repmat(model.MUAprdRest(2, :).*0.2, nSynapses, 1)';

%% Step 2: generate dendritic current

synapseFunc = @(x) zeromean(2*rand(x,1)-1); 
%synapseFunc = @(x) rand(x, 1);
totalSpikes = [];
nTrials = 150;

for k = 1 : 2
    for k1 = 1 : nTrials
        tmp    = rand(size(rate{k})); 
        spikes = zeros(size(tmp));
        
        spikes(tmp < rate{k}) = 1;
        
        peakCurrent = synapseFunc(nSynapses);
        spikes      = bsxfun(@times, spikes, peakCurrent');
        
        % sum over synaoses
        totalSpikes{k}(k1, :) = sum(spikes, 2);
    end
end

figure (3), clf
subplot(2, 2, 1)
plot(model.MUAprdRest')

subplot(2, 2, 2)
plot(abs(sum(totalSpikes{1})))

%% Compute post synaptic and dendritic current

% two time constants:
% time constant for dendritic integration
alpha = 0.1; % from Miller et al.

% time constant for post-synaptic current
tau = 0.0023; % from Miller et al.

Q = [];
I = {};
I{1} = zeros(nTrials, nt);
I{2} = zeros(nTrials, nt);

psc = exp(-1/tau*(0:dt:.100));

% convolve spikes with post synaptic current
for k = 1 : 2
    for k1 = 1 : nTrials
        Q{k}(k1, :) = convCut(totalSpikes{k}(k1, :), psc, nt);
    end
end

m = zeros(1, length(Q));
for k = 1 : 2
    for k1 = 1 : nTrials
        for jj = 1 : length(Q{1}) - 1
            % rate of change in current
            dIdt = (Q{k}(k1, jj) - I{k}(k1, jj)) / alpha;
            
            % stepwise change in current
            dI = dIdt * dt;
            
            % current at next time point
            I{k}(k1, jj+1) = I{k}(k1, jj) + dI;
        end
    end
end

%% Extract broadband

extractedbb = {};

for k = 1 : 2
    for k1 = 1 : nTrials
        extractedbb{k}(k1, :) = extractBroadband(I{k}(k1, :), srate, 4, bands);
    end
end

%%
% figure (4), clf
% subplot(2, 2, 1), plot(Q{1}'), title('Q')
% subplot(2, 2, 2), plot(I{1}'), title('I')
% subplot(2, 2, 4)
% plot(t(crange), mean(extractedbb{1}))

figure (2), 
subplot(2, 2, 1), plot(t(crange), normMax(median(extractedbb{1})), 'r-')
subplot(2, 2, 2), plot(t(crange), normMax(median(extractedbb{2})), 'r-')