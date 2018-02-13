function totalSpikes = dn_generatePoissonSpikes(spikeRate, nTrials, nSynapse)

% adapted from Jon Winawer's 2013 Current Biology paper

%% useful functions

synapseFunc = @(x) zeromean(2 * rand(x,1)-1);

%%

for k = 1 : nTrials
    tmp = rand(size(spikeRate)); spikes = zeros(size(tmp));
    
    spikes(tmp < spikeRate) = 1;
    
    peakCurrent = synapseFunc(nSynapse);   % multiply each spike by the peak current of the appropriate synapse
    
    spikes      = bsxfun(@times, spikes, peakCurrent');
    
    totalSpikes(k, :) = sum(spikes, 2); % output spikes
end


end