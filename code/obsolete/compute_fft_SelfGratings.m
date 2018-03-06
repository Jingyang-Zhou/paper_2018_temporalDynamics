% inspect the power spectrum

%% load data

dataLoc = '/Volumes/server/Projects/Temporal_integration/data/ECoG/';
dataNm1 = 'CONTEXT.mat';

normMax = @(x) (x - mean(x(1 : 100)))./max(x - mean(x(1 : 100)));
normMax_noBase = @(x) x./max(x);

a = load(fullfile(dataLoc, dataNm1));

mua = a.MUA;
lfp = a.LFP;

m_mua = median(mua);
m_lfp = median(lfp);

dt = 1/a.FsD;
t  = (0:size(lfp, 2)-1)*dt;

figure (1), clf
subplot(2, 2, 1)
plot(t, m_mua, 'k-'), hold on, ylim([2, 2.6])


%% Compute power spectrum of the LFP signal

f_lfp  = [];
f_lfp1 = [];

epochLength = .4;
baselineIdx = t<epochLength;
stimIdx     = t > .45 & t -.45-dt < epochLength;
t1 = t(baselineIdx);
t2 = t(stimIdx);
t2 = t2 - t2(1);
nf = sum(stimIdx);

for k = 1 : size(lfp, 1)
    f_lfp(k, :)  = abs(fft(lfp(k, baselineIdx)));
    f_lfp1(k, :) = abs(fft(lfp(k, stimIdx)));
end

f  = [0  : nf - 1]./ max(epochLength);

fm_lfp = geomean(f_lfp);
fm_lfp1 = geomean(f_lfp1);

%%
subplot(2, 2, 2)
loglog(f, fm_lfp, 'k-'), hold on,
loglog(f, fm_lfp1, 'r-')
axis tight, box off
set(gca, 'fontSize', 14), legend('no stim', 'stim on')
xlim([0, 200]);

subplot(2, 2, 3)
plot(f, fm_lfp1./fm_lfp), xlim([0, 200]);

%% compute spectrogram
spectro = [];
spectrbase = [];
params.pad=-1;
params.tapers=[3 5];
params.fpass=[0 200];
params.Fs= a.FsD;
params.trialave=1;

movingwin=[.100 .01];

for k = 1 : size(lfp, 1)
    [spectrobase(k, :,:), t_tf, f] = mtspecgramc(lfp(k, baselineIdx), movingwin, params);
    [spectro(k, :,:), t_tf, f] = mtspecgramc(lfp(k, stimIdx), movingwin, params);
end

spectro = spectro./spectrobase;
m_spectro = squeeze(mean(spectro));

subplot(2, 2, 4)
imagesc(t_tf-0.05, f, log10(m_spectro)')
axis xy, set(gca, 'fontSize', 14),
xlabel('time (s)'), ylabel('frequency (Hz)')

%% check Hermes 2015 Cerebral Cortex papers

gamOscPth = '/Users/winawerlab/matlab/git/Papers_Hermes_2015_CerCortex';
addpath(genpath(gamOscPth))

