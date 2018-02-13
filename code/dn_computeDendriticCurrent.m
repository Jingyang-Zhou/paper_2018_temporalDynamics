function I = dn_computeDendriticCurrent(totalSpikes, tau, alpha, dt)

% This function is based on ecogSimulate.m from Jon Winawer's 2012 Current
% Biology paper. 

% INPUTS: 

% totalSpikes: of dimension [number of trials x time course].
% tau        : time constant for post-synaptic current
% alpha      : time constant for dendritic integration

% OUTPUTS: 

% I : of dimension [number of trials x time course], output current
% voltage?

%% prep

nTrials = size(totalSpikes, 1);
nt      = size(totalSpikes, 2);

%% generate post-synaptic and dendritic current

Q = [];
I = zeros(nTrials, nt);

% compute post-synaptic current
psc = exp(-1/tau*(0:dt:.100));

for k = 1 : nTrials
    Q(k, :) = convCut(totalSpikes(k, :), psc, nt);
    for jj = 1 : nt - 1
        % rate of change in current
        dIdt = (Q(k, jj) - I(k, jj))./alpha;
        
        % stepwise change in current
        dI = dIdt * dt;
        
        % current at next time point
        I(k, jj+1) = I(k, jj) + dI;
    end
end


end