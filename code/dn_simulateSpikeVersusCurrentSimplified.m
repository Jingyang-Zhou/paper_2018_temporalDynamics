% QUESTION:

% The time constants estimated in neuronal responses are very
% different for single unit and ECoG measurements. The hypothesis is that
% the difference could be caused by that ECoG is measuring the dentritic
% currents and single units are measuring neuronal spiking activity. This
% simulation is to verify this idea.

% Currently I am assuming one population of synapses instead of 2
% (as in Jon's previous simulation for the 2013 current biology paper.)

temporalPth = '/Volumes/server/Projects/Temporal_integration/code/';
addpath(genpath(temporalPth))

simPth      = '~/Google Drive/data/ECoGBOLD2013/';
addpath(genpath(simPth))

bbPth       = '~/Google Drive/broadband tutorial';
addpath(genpath(bbPth))

%% useful functions

normMax = @(x) x./max(x);
normMax_baseline = @(x) (x - mean(x(1 : 150)))./max(x - mean(x(1 : 150)));

%% pre-defined variables and single unit repsonses/spike rate

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
shift  = 0.05;  % time depay between stimulus onset and response onset
scale  = 1;     % max response amplitude?

% resting state spike rate:
restingRate = 0.1;

% compute the normalization model (approximates single unit measurement)
rate = normMax(dn_DNmodel([tau1, weight, tau2, n, sigma, shift, scale], stim, t));
rate = normMax((rate + restingRate))./5;
%
% two time constants:
% time constant for dendritic integration
alpha = 0.1; % from Miller et al.

% time constant for post-synaptic current
tau = 0.0023; % from Miller et al.

nTrials = 20;

%% Derived variables

nt = length(t) * nTrials; % number of samples
dt = median(diff(t)); % step size in ms

%% simulate dendritic current responses (fast way)

ts = [];

%  Stage (1) Extend Poisson spike rates:

% The rate needs to be a matrix num_time_points x num_synapses. If the rate
% is a scalar, it needs to be expanded in two dimensions. If the rate is a
% vector (ie time-varying), it needs to be expanded only for the number of
% synapses.
if numel(rate) == 1,  n = nt; else n = 1; end
rate  = repmat(rate(:),  n, numSynapses);

synapseFunc = @(x) zeromean(2*rand(x,1)-1); % uniformly distributed random variable centered around 0
%synapseFunc = @(x) rand(1, x);


%% generate dendritic current

totalSpikes = [];

% generate spike rates:
for k = 1 : nTrials
    spikes = [];
    % % We could use a true Poisson process (slow)
    % spikes = poissrnd(rate)... or an approximation (faster)
    tmp = rand(size(rate)); spikes = zeros(size(tmp));
    spikes(tmp < rate) = 1;
    
    % multiply each spike by the peak current of the appropriate synapse
    peakCurrent = synapseFunc(numSynapses);
    spikes = bsxfun(@times, spikes, peakCurrent);
    
    % sum over synapses
    totalSpikes(:, k) = sum(spikes,2);
end

totalSpikes = reshape(totalSpikes, [numSynapses * nTrials, 1]);

%% Compute post synaptic and dendritic current

I = zeros(1, nt);
Q = [];

psc = exp(-1/tau*(0:dt:.100));

% convolve spikes with post synaptic current
Q = conv(totalSpikes, psc, 'full'); Q = Q(1:nt);

for jj = 1:length(Q)-1
    
    % rate of change in current
    dIdt = (Q(jj) - I(jj)) / alpha;
    
    % stepwise change in current
    dI = dIdt * dt;
    
    % current at next time point
    I(jj+1) = I(jj) + dI;
end

%% extract the broadband signal

bb = [];
output = [];
srate = 1000;

% band_rg  = [80 200];
% band_w   = 20;
% lb       = band_rg(1):band_w:band_rg(2)-band_w;
% ub       = lb+band_w;
% bands   = [lb; ub]';

bands = {[70, 170], 20};


for k = 1 : nTrials
    bb = extractBroadband(I, srate, 4, bands);
   % bb = extract
end

bb = reshape(bb, [numSynapses, nTrials]);
rsI = reshape(I, [numSynapses, nTrials]);

%% normalize the two time series

mbb   = normMax_baseline(median(bb, 2));
mrate = normMax_baseline(median(rate, 2));
mrsI  = median(rsI, 2);

%% visualize

figure (1), clf
subplot(2, 2, 1), plot(t, rate, 'k-'), hold on, plot(t, stim.*max(rate), 'k:'), box off
xlabel('time (s)'), ylabel('poisson mean'), title('poisson mean for spike rates')

subplot(2, 2, 2), stem(t, abs(totalSpikes(1 : length(t))), 'k-',  'markersize', 1), box off, xlabel('time'), title('spike rates per trial')

subplot(2, 2, 3), plot(t, mrsI), title('median time course in one trial'),  box off, xlabel('time (s)'), ylabel('voltage')

subplot(2, 2, 4), plot(t, mbb, 'k-', 'linewidth', 2), hold on, plot(t, mrate, 'r-', 'linewidth', 2),
legend('dendritic current', 'spike rate'), box off, xlabel('time'), ylabel('amplitude')
