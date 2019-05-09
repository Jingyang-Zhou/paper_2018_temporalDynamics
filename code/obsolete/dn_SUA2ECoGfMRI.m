% dn_from dingle unit to ECoG and fMRI

%% useful functions

normMax_range = @(x, range) (x - mean(x(range)))./max(x - mean(x(range)));
normMax       = @(x) x./max(x);

%% load single unit data

dataLoc = fullfile(temporalRootPath, 'data', 'ECoG');
fName   = 'figure4Data.xlsx';

a = xlsread(fullfile(dataLoc, fName));

t      = a(:, 1)./1000;
data   = a(:, 2 : end)';
mdata  = normMax_range(median(data, 1), 1 : 20);

nCells = size(data, 1);

cellNum = [1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 3];

whichcells = 1 : 12;

%% load ecog data

dataLoc = fullfile(temporalRootPath, 'data');
fName   = 'dn_data.mat';

a       = load(fullfile(dataLoc, fName));
param   = a.param;

%% extend data and prep for bootstraps
extend_data = data(1, :);

for k = 1 : nCells
    tmp = repmat(data(k, :), [cellNum(k), 1]);
    extend_data = [extend_data; tmp];
end

%% bootstrap:
nBoots    = 100;
boot_data = [];

for k = 1 : nBoots
    idx = randi(sum(cellNum), [1, sum(cellNum)]);
    boot_data(k, :) = normMax(mean(extend_data(idx, :)));
end

%% model grid fit to single unit data

stim = ones(1, length(t));
init = [0.08, 0.08, 2, 0.1, 0];
lb   = [0.01,0.01,0,0.01,0];
ub   = [1, 1, 6, 0.5, 0.1];
sua.param = [];

sua = [];

% for k = 1 : nBoots
%     k
%     [sua.seed(k, :), sua.seedR2(k)] = dn_gridFit(boot_data(k, :), param, stim, t, 'uniphasic');
% end

sua.seed = repmat([0.05, 0.05, 2, 0.01, 0], [nBoots, 1]);
%% model fine fit to single unit data

for k = 1 : nBoots
    [sua.param(k, :), ~, exitflg(k)] = fminsearchbnd(@(x) dn_computeFineFit(x, boot_data(k, :), stim, t, 'uniphasic'), init, lb, ub);
    sua.pred(k, :)  = normMax(dn_DNmodel([sua.param(k, 1), 0, sua.param(k, 2 : end), 1], stim, t));
end

% compute summary metrics:
sua.summary = dn_computeDerivedParams(sua.param, 'uniphasic');


%% plot single unit data and model prediction

sua.msummaryPrm = [median(sua.summary.t2pk), median(sua.summary.r_asymp)];
sua.ssummaryPrm = [prctile(sua.summary.t2pk, [25, 75]); prctile(sua.summary.r_asymp, [25, 75])];

eocg.msummaryPrm =[median(a.dn.derivedPrm.t2pk(1 : 100)), median(a.dn.derivedPrm.r_asymp(1 : 100))];
eocg.msummaryPrm =[prctile(a.dn.derivedPrm.t2pk(1 : 100), [25, 75]); prctile(a.dn.derivedPrm.r_asymp(1 : 100), [25, 75])];

figure (1), clf

for k = 1 : length(whichcells)
    subplot(4, 3, k),
    plot(t, data(whichcells(k), :), 'k-', 'linewidth', 3), hold on, % plot(t, sua.pred(k, :), 'r-')
    box off, title(cellNum(k)), set(gca, 'xtick', [0, 0.2], 'ytick', [0, 1], 'fontsize', 12), axis square
end

figure (2), clf
yyaxis left,plot(1, sua.msummaryPrm(1), 'ro', 'markersize', 10, 'markerfacecolor', 'r'), hold on,
plot(1, eocg.msummaryPrm(1), 'ko', 'markersize', 10, 'markerfacecolor', 'k'), ylim([0, 300])
plot([1, 1], eocg.msummaryPrm(1, :), 'k-', 'linewidth', 2)
plot([1, 1], sua.ssummaryPrm(1, :), 'k-', 'linewidth', 2), ylabel('t2pk (ms)')

yyaxis right,plot(2, sua.msummaryPrm(2), 'ro',  'markersize', 10, 'markerfacecolor', 'r'),
plot(2, eocg.msummaryPrm(2), 'ko', 'markersize', 10, 'markerfacecolor', 'k'), %ylim([0, 0.5]),
plot([2, 2], sua.ssummaryPrm(2, :), 'k-', 'linewidth', 2)
plot([2, 2], eocg.msummaryPrm(2, :), 'k-', 'linewidth', 2), ylabel('r-asymp')
set(gca, 'xtick', [1, 2], 'xticklabel', {'t2pk', 'r-asymp'}, 'fontsize', 16), box off
xlim([0.5, 2.5])

figure (3), clf % plot bootstrapped data and model fit
boundedline(t, mean(boot_data), std(boot_data), 'k-'), hold on
%boundedline(t, mean(sua.pred), std(sua.pred), 'r-')
plot(t, median(sua.pred), 'r-', 'linewidth', 3), axis square, axis tight, set(gca, 'xtick', [0, 0.2], 'ytick', [0, 1], 'fontsize', 18)
xlabel('time (s)'), ylabel('normalized response')

%% average parameters

lfp = [];

t1    = 0.001 : 0.001 : 0.9;
stim1 = zeros(1, length(t1));
stim1(1 : 200) = 1;

lfp.param = median(sua.param);

lfp.spikeRate = normMax(dn_DNmodel([lfp.param(1), 0, lfp.param(2 : end), 1], stim1, t1));

%% FROM MUA TO LFP

spikeRate = [];

restSpikeRate = 0.05;

lfp.spikeRateRest = normMax(lfp.spikeRate + restSpikeRate);

nSynapses = 1000;

rate = repmat(lfp.spikeRateRest.*0.05, nSynapses, 1)';

%% generate dendritic current

synapseFunc = @(x) zeromean(2*rand(x,1)-1);
totalSpikes = [];
nTrials     = 100;

for k = 1 : nTrials
    tmp    = rand(size(rate));
    spikes = zeros(size(tmp));
    
    spikes(tmp < rate) = 1;
    
    peakCurrent = synapseFunc(nSynapses);
    spikes      = bsxfun(@times, spikes, peakCurrent');
    
    % sum over synapses
    totalSpikes(k, :) = sum(spikes, 2);
end

figure (4), clf
plot(mean(totalSpikes))

%% compute post-synaptic and dendritic current
% two time constants:
% time constant for dendritic integration
alpha = 0.1; % from Miller et al.

% time constant for post-synaptic current
tau = 0.0023; % from Miller et al.

dt = 0.001;

Q = [];
I = zeros(nTrials, length(t1));

psc = exp(-1/tau*(0:dt:.100));

for k = 1 : nTrials
    Q(k, :) = convCut(totalSpikes(k, :), psc, length(t1));
end

for k = 1 : nTrials
    for jj = 1 : length(Q) - 1
        % rate of change in current
        dIdt = (Q(k, jj) - I(k, jj)) / alpha;
        
        % stepwise change in current
        dI = dIdt * dt;
        
        % current at next time point
        I(k, jj+1) = I(k, jj) + dI;
    end
end

%% extract broadband
srate = 1000;
bands = {[70 170], 20};

bb = [];

bb = extractBroadband(I', srate, 4, bands);

mbb = normMax(mean(bb, 2));

figure (5), clf
plot(t1, lfp.spikeRate, 'r-', 'linewidth', 3), hold on
plot(t1, mbb, 'b-', 'linewidth', 3), xlim([0, 0.2]), axis square, box off
set(gca, 'xtick', [0, 0.2], 'ytick', [0, 1], 'fontsize', 18)
xlabel('time (s)'), ylabel('normalized repsonse')
