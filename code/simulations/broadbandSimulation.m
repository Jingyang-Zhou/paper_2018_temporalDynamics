% broadband simulation

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

%% PRE-DEFINED VARIABELS AND FUNCTIONS

% USEFUL FUNCTIONS --------------------------------------------------------
whiten   = @(x) (x - mean(x(:)))./ diff(prctile(x, [.25 .75]));
normMax  = @(x) x./max(x);
normBase = @(x) x-mean(x(1 : 200));

% PREDEFINED VARIABLES ----------------------------------------------------
srate    = 1000;
nTrials  = 50;

% HILBERT TRANSFORMS ------------------------------------------------------
stand_hilbert = @(x) abs(hilbert(x));

% DN MODEL PARAMETERS -----------------------------------------------------
dnParams = [0.1, 0, 0.1, 2, 0.1, 0, 1];
t        = 0.001 : 0.001 : 1.2;
stim     = zeros(1, length(t));
stim(t>=0.2 & t<0.7) = 1;

%% FORWARD MODEL (BUILD THE GROUND TRUTH)
vol = [];

for k = 1 : nTrials
    [vol(k, :), dn] = dn_LFP_forwardModel(dnParams, t, stim);
end

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

bands{5} = {[20, 200], 20};
bands{6} = {[100, 200], 20};

% EXTRACT AND COMPARE BROADBAND USING THE ABOVE BANDWIDTH ----------------
bp = {}; bb = {}; m_bb = [];

for k = 1 : length(bands)
    nBands     = (bands{k}{1}(2) - bands{k}{1}(1))/bands{k}{2};
    bp{k}      = squeeze(bandpassFilter(vol, srate, bands{k}));
    % take Hilbert transform 
    bb{k}      = squeeze(geomean(stand_hilbert(bp{k}).^2, 3));
    m_bb(k, :) = normMax(mean(bb{k}, 1));
end

% VISUALIZE --------------------------------------------------------------
figure 
for k = 1 : length(bands)
    subplot(6, 1, k)
    plot(t, m_bb(k, :), 'r-', 'linewidth', 2), hold on
    plot(t, normMax(dn), 'k-', 'linewidth', 2), box off
end

