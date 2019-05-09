% parameter sweeping

% tau1  : 0.07 : 1
% tau2  : 0.07 : 1
% n     : 1 : 6
% sigma : 0.01 : 0.5
% w     : 0 : 1

%% MAKE STIMULUS

t    = 0.001 : 0.001 : 1.2;
stim = zeros(size(t));
stim(201 : 700) = 1;

stim_low = stim .*0.3;

normMax = @(x) x ./ max(x);

%% PARAMETER SWEEPING RESULTS

rsp    = [];
prm    = [0.07, 0, 0.1, 2, 0.1, 0, 1];
nsteps = 10;

tau1  = round(linspace(0.07, 1, nsteps), 2);
tau2  = round(linspace(0.07, 1, nsteps), 2);
n     = round(linspace(1, 6, nsteps), 2);
sigma = round(linspace(0.01, 0.5, nsteps), 2);
w     = round(linspace(0, 1, nsteps), 2);

% changing tau1
for k = 1 : nsteps
    prm1    = prm;
    prm1(1) = tau1(k);
    rsp(1, k, 1, :) = dn_DNmodel(prm1, stim, t);
    rsp(1, k, 2, :) = dn_DNmodel(prm1, stim_low, t);
    tmp = max(squeeze(rsp(1, k, :)));
    rsp(1, k, :, :) = rsp(1, k, :, :)./tmp;
end

% changing tau2
for k = 1 : nsteps
    prm2    = prm;
    prm2(3) = tau2(k);
    rsp(2, k, 1, :) = dn_DNmodel(prm2, stim, t);
    rsp(2, k, 2, :) = dn_DNmodel(prm2, stim_low, t);
    tmp = max(squeeze(rsp(2, k, :)));
    rsp(2, k, :, :) = rsp(2, k, :, :)./tmp;
end

% changing n
for k = 1 : nsteps
    prm3    = prm;
    prm3(4) = n(k);
    rsp(3, k, 1, :) = dn_DNmodel(prm3, stim, t);
    rsp(3, k, 2, :) = dn_DNmodel(prm3, stim_low, t);
    tmp = max(squeeze(rsp(3, k, :)));
    rsp(3, k, :, :) = rsp(3, k, :, :)./tmp;
end

% changing sigma
for k = 1 : nsteps
    prm4    = prm;
    prm4(5) = sigma(k);
    rsp(4, k, 1, :) = dn_DNmodel(prm4, stim, t);
    rsp(4, k, 2, :) = dn_DNmodel(prm4, stim_low, t);
    tmp = max(squeeze(rsp(4, k, :)));
    rsp(4, k, :, :) = rsp(4, k, :, :)./tmp;
end

% changing w
for k = 1 : nsteps
    prm5    = prm;
    prm5(2) = w(k);
    rsp(5, k, 1, :) = dn_DNmodel(prm5, stim, t);
    rsp(5, k, 2, :) = dn_DNmodel(prm5, stim_low, t);
    tmp = max(squeeze(rsp(5, k, :)));
    rsp(5, k, :, :) = rsp(5, k, :, :)./tmp;
end

%% PARAMETER SWEEPING FIGURE

figure (1), clf
for k = 1 : nsteps
    subplot(5, nsteps, k)
    plot(t, squeeze(rsp(1, k, :, :)), 'k-'), hold on, axis tight, ylim([-0.1, max(rsp(:))]), hold on, box off, set(gca, 'ytick', '', 'xtick', ''), title(tau1(k))
    subplot(5, nsteps, k + nsteps)
    plot(t, squeeze(rsp(2, k, :, :)), 'k-'), hold on, axis tight, ylim([-0.1, max(rsp(:))]), hold on, box off, set(gca, 'ytick', '', 'xtick', ''), title(tau2(k))
    subplot(5, nsteps, k + nsteps * 2)
    plot(t, squeeze(rsp(3, k,:, :)), 'k-'), hold on, axis tight, ylim([-0.1, max(rsp(:))]), hold on, box off, set(gca, 'ytick', '', 'xtick', ''), title(n(k))
    subplot(5, nsteps, k + nsteps * 3)
    plot(t, squeeze(rsp(4, k, :,:)), 'k-'), hold on, axis tight, ylim([-0.1, max(rsp(:))]), hold on, box off, set(gca, 'ytick', '', 'xtick', ''), title(sigma(k))
    subplot(5, nsteps, k + nsteps * 4)
    plot(t, squeeze(rsp(5, k,:, :)), 'k-'), hold on, axis tight, ylim([-0.1, max(rsp(:))]), hold on, box off, set(gca, 'ytick', '', 'xtick', [0.2, 0.7], 'xticklabel', [0, 0.5]), title(w(k))
end

for k = 1 : nsteps * 5
    subplot(5, nsteps, k), plot([0.2, 0.2], [0, 1], 'r:'), plot([0.7, 0.7], [0, 1], 'r:')
end
