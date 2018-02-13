function [I, dn] =  dn_LFP_forwardModel(dnParams, t, stim)

% INPUTS --------------------------------------------------------------
% dnParams: DN model parameters
% dnStim  : stimulus input to the DN model
% t       : time 

% OUTPUTS -------------------------------------------------------------

%% EXAMPLES 

dnParams = [0.1, 0, 0.1, 2, 0.1, 0, 1];
dt       = 0.001;
t        = dt : dt : 1.2;
stim     = zeros(1, length(t));
stim(t>=0.2 & t<0.7) = 1;

figureOn = 1;

%% PRED-DEFINED PARAMETERS AND FUNCTIONS

tau         = 0.0023; % time scale parameter for post-synaptic current
alpha       = 0.05; % time constant for dendritic integration 
srate       = 1000;
restingRate = 0.1;
%nTrials     = 100;
nSynapses   = 1000;
 dt         = 0.001;

normMax     = @(x) x./max(x);
synapseFunc = @(x) zeromean(2*rand(x,1)-1); % uniformly distributed random variable centered around 0

% DERIVED PARAMETERS ----------------------------
t_lth   = length(t);

%% GENERATE DN MODEL PREDICTION AND SPIKE RATES

dn     = normMax(normMax(dn_DNmodel(dnParams, stim, t)) + restingRate);
weight = synapseFunc(nSynapses);

totalSpikes = [];

for k = 1 %: nTrials
    tmp = rand(size(dn)); spikes = zeros(size(tmp));
    spikes(tmp < dn) = 1;
    
    % multiply each spike by the peak current of the appropriate synapse
    peakCurrent = synapseFunc(nSynapses);
    spikes      = bsxfun(@times, spikes, peakCurrent);
    % sum over synapses
    totalSpikes(:, k) = sum(spikes);
end

%% LEAKY INTEGRATION

psc = exp(-1/tau*(0:dt:.100));

I = zeros(1, t_lth);
Q = [];

psc = exp(-1/tau*(0:dt:.100));

% convolve spikes with post synaptic current
Q = conv(totalSpikes, psc, 'full'); Q = Q(1:t_lth);

for jj = 1:length(Q)-1
    
    % rate of change in current
    dIdt = (Q(jj) - I(jj)) / alpha;
    
    % stepwise change in current
    dI = dIdt * dt;
    
    % current at next time point
    I(jj+1) = I(jj) + dI;
end

%% VISUALIZE

if figureOn == 1
    figure (100), clf
    plot(t, I), hold on
    
    %plot(t, dn, 'k:')
end

end