% dn_mkFigure4_illustration

%% make 3 time course prediction using the biphasic DN model

%% set parameters

normMax = @(x) x./max(x);

weight = [0, 0.6, 0.75];

stim = zeros(1, 1200);
stim(200 : 700) = 1;

t = [1 : length(stim)]./1000;

% fields = {'tau1', 'weight', 'tau2', 'n', 'sigma', 'shift', 'scale'};

for k = 1 : length(weight)
    param = [0.05, weight(k), 0.057, 2, 0.1, 0, 1];
    normrsp(k, :) = normMax(dn_DNmodel(param, stim, t));
end

%% plot 

figure (4), clf
patch([0.2, 0.7, 0.7, 0.2], [0, 0, 1, 1], 0.8 * ones(1, 3)), hold on
plot(t, normrsp', 'k-', 'linewidth', 3), box off
set(gca, 'xtick', [0.2, 0.7], 'xticklabel', [0, 0.5])
set(gca, 'ytick', [0, 1], 'fontsize', 14)