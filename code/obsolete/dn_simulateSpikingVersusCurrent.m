% QUESTION:

% The time constants estimated in neuronal responses are very
% different for single unit and ECoG measurements. On hypothesis is that
% the difference could be caused by ECoG is measuring the dentritic
% currents and single units are measuring neuronal spiking activity. This
% simulation is to verify this idea.

% Currently I am assuming one population of synapses instead of 2
% (as in Jon's previous simulation for the 2013 current biology paper.)

temporalPth = '/Volumes/server/Projects/Temporal_integration/code/';
addpath(genpath(temporalPth))

simPth      = 'Google Drive/data/ECoGBOLD2013/';
addpath(genpath(simPth))

%% useful functions

normMax = @(x) x./max(x);
normMax_baseline = @(x) (x - mean(x(1 : 150)))./max(x - mean(x(1 : 150)));

%% pre-defined variables and single unit repsonss/spike rate

rate = []; % poisson spike rate
t    = 0.001 : 0.001 : 1;

% make stimulus
stim = zeros(1, length(t));
stim(t>0.2 & t<=0.7) = 1;

% pre-define number of synapses
numSynapses = 1000;

% Poisson rate for the population:

% normalization model parameters
tau1   = 0.007; % summation time constant
weight = 0;     % weight of the second pulse in the biphasic IRF
tau2   = 0.01;  % the adaptation time constant
n      = 2;     % normalization power
sigma  = 0.05;  % normalization semi-saturation
shift  = 0.05;     % time depay between stimulus onset and response onset
scale  = 1;     % max response amplitude?

% resting state spike rate:
restingRate = 0.1;

% compute the normalization model (approximates single unit measurement)
rate = normMax(dn_DNmodel([tau1, weight, tau2, n, sigma, shift, scale], stim, t));
rate = normMax((rate + restingRate));

% two time constants:
% time constant for dendritic integration
alpha = 0.1; % from Miller et al.

% time constant for post-synaptic current
tau = 0.0023; % from Miller et al.

%% Derived variables

nt = length(t); % number of samples
dt = median(diff(t)); % step size in ms

%% visualize single unit spike rate measurement

figure (1), clf
subplot(2, 2, 1), plot(t, rate, 'k-'),hold on, plot(t, stim, 'k:'),  box off
xlabel('time (s)'), ylabel('amplitude'), title('Poisson spike rate over time'), legend('spike rate', 'stimulus')

%% COMPUTE DENDRITIC CURRENT:

%%  Stage (1) Extend Poisson spike rates:

% The rate needs to be a matrix num_time_points x num_synapses. If the rate
% is a scalar, it needs to be expanded in two dimensions. If the rate is a
% vector (ie time-varying), it needs to be expanded only for the number of
% synapses.
if numel(rate) == 1,  n = nt; else n = 1; end
rate  = repmat(rate(:),  n, numSynapses);

%% Stage (2) Generate spike rate

spikes = [];

% % We could use a true Poisson process (slow)
% spikes = poissrnd(rate)... or an approximation (faster)
tmp = rand(size(rate)); spikes = zeros(size(tmp));
spikes(tmp < rate) = 1;

%% Stage (3) Generate currents

synapseFunc = @(x) zeromean(2*rand(x,1)-1);
%synapseFunc = @(x) rand(1, x);
peakCurrent = synapseFunc(numSynapses);

% multiply each spike by the peak current of the appropriate synapse
spikes = bsxfun(@times, spikes, peakCurrent);

% sum over synapses
spikes = sum(spikes,2);

% Post-synaptic current. This current (multiplied by the synaptic weight)
% is initiated each time a spike arrives. It rises quickly and falls
% slowly.
psc = exp(-1/tau*(0:dt:.100));

% convolve spikes with post synaptic current
Q = conv(spikes, psc, 'full'); Q = Q(1:nt);

figure (1), subplot(2, 2, 2)
stem(t, spikes(:, 1), 'k-', 'markersize', 1), hold on,
plot(t, stim.*200, 'r:', 'linewidth', 3), box off, xlabel('time (s)')
title('spike rates')

%% Leaky integration

% Initialize the dendritic current
I = zeros(1, nt);

% We calculate the stepwise change in current because the decay is
% proportional to the current level (described by a differential
% equation)
for jj = 1:length(Q)-1
    
    % rate of change in current
    dIdt = (Q(jj) - I(jj)) / alpha;
    
    % stepwise change in current
    dI = dIdt * dt;
    
    % current at next time point
    I(jj+1) = I(jj) + dI;
end

% keep only time points with t > 0. time points less than zero indicate a
% 'pre-conditioning' that is presumably used to let the system settle.

if min(t) < 0,
    inds = t>0;
    t = t(inds);
    I = I(inds);
end

% If no outputs, then plot
%figure(1); subplot(2, 2, 4), plot(t, I);

%% extract the broadband signal

F = fft(I);

band_rg  = [80 200];
band_w   = 10;
lb       = band_rg(1):band_w:band_rg(2)-band_w;
ub       = lb+band_w;
bands   = [lb; ub]';

%    filter
whiten = @(x) bsxfun(@plus, zscore(x), mean(x));

srate = 1000;

signal.bp_multi = zeros(length(I), size(bands, 1));
for ii = 1 : size(bands, 1)
    [signal.bp_multi(:,ii)]=butterpass_eeglabdata(I',bands(ii,:), srate);
end


signal.bb = normMax_baseline(sum(whiten(abs(hilbert(signal.bp_multi))), 2));
norm_rate = normMax_baseline(mean(rate, 2));

figure (1),
subplot(2, 2, 3), plot(t, I),
subplot(2, 2, 4), plot(t, signal.bb, 'k-', 'linewidth', 2), hold on, plot(t, norm_rate, 'r-', 'linewidth', 2),
legend('ECoG broadband', 'spike rate'), box off, xlabel('time'), ylabel('amplitude')
