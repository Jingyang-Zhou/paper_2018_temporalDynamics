% simualate_ static versus dynamics normalization predictions to time to
% peak

%% MAKE CONTRAST STIMULI

t = 0.001 : 0.001 : 1.5;

nContrast = 10;
contrastLevels = linspace(0.1, 1, nContrast);

nStim = 10;
stim  = zeros(nStim, length(t));

for k = 1 : nContrast
   stim(k, 200 : 700) = contrastLevels(k); 
end

% VISUALIZE STIMULI
figure (1), clf, set(gca, 'colororder', copper(nContrast)), hold on
plot(t, stim, 'linewidth', 2), 
xlabel('time (s)'), title('stimuli')

%% MAKE LINEAR MODEL PREDICTION

tau1 = 0.1;
irf = gammaPDF(t, tau1, 2);

% MAKE LINEAR REPSONSE TO ALL STIMULI

for k = 1 : nContrast
    prd.lin(k, :) = convCut(irf, stim(k, :), length(t));
    prd.linMaxIdx(k) = find(prd.lin(k, :) == max(prd.lin(k, :)));
end

figure (2), 
subplot(2, 3, 1),cla, set(gca, 'colororder', copper(nContrast)), hold on
plot(t, prd.lin), axis tight
for k = 1 : nContrast
    plot(t(prd.linMaxIdx(k)), prd.lin(k, prd.linMaxIdx(k)), 'ro', 'markerfacecolor', 'r')
end
title('Linear model'), set(gca, 'fontsize', 14)

subplot(2, 3, 4), cla
plot(contrastLevels, t(prd.linMaxIdx), 'ro:', 'markerfacecolor', 'r'), hold on
plot(contrastLevels, prd.lin(:, prd.linMaxIdx(k)), 'bo:', 'markerfacecolor', 'b'), box off

%% MAKE CTS MODEL PREDICTION

n     = 2;
sigma = 0.1;

cts = @(x, n, sigma) x.^n./(sigma^n + x.^n); % x is the linear response

for k = 1 : nContrast
   prd.cts(k, :)    = cts(prd.lin(k, :), n, sigma);
   prd.ctsMaxIdx(k) = find(prd.cts(k, :) == max(prd.cts(k, :)));
end

figure (2), subplot(2, 3, 2), cla, set(gca, 'colororder', copper(nContrast)), hold on
plot(t, prd.cts), hold on
for k = 1 : nContrast
   plot(t(prd.ctsMaxIdx(k)), prd.cts(k, prd.ctsMaxIdx(k)), 'ro', 'markerfacecolor', 'r') 
end
title('CTS model'), set(gca, 'fontsize', 14)

subplot(2, 3, 5), cla, 
plot(contrastLevels, t(prd.ctsMaxIdx), 'ro:', 'markerfacecolor', 'r'), hold on
plot(contrastLevels, prd.cts(:, prd.ctsMaxIdx(k)), 'bo:', 'markerfacecolor', 'b'), box off

%% MAKE DN MODEL PREDICTIONS
tau2   = 0.1;
dnPrms = [tau1, 0, tau2, n, sigma, 0, 1];

for k = 1 : nContrast
   prd.dn(k, :) = dn_DNmodel(dnPrms, stim(k, :), t); 
   prd.dnMaxIdx(k) = find(prd.dn(k, :) == max(prd.dn(k, :)));
end

prd.dn = prd.dn ./max(prd.dn(:));

figure (2), subplot(2, 3, 3), cla, set(gca, 'colororder', copper(nContrast)), hold on
plot(t, prd.dn), hold on
for k = 1 : nContrast
   plot(t(prd.dnMaxIdx(k)), prd.dn(k, prd.dnMaxIdx(k)), 'ro', 'markerfacecolor', 'r')
end
title('DN model'), set(gca, 'fontsize', 14)

subplot(2, 3, 6), cla, 
plot(contrastLevels, t(prd.dnMaxIdx), 'ro:', 'markerfacecolor', 'r'), hold on
plot(contrastLevels, prd.cts(:, prd.dnMaxIdx(k)), 'bo:', 'markerfacecolor', 'b'), box off, 
legend('time to peak', 'peak amplitude', 'Location', 'best'), xlabel('contrast levels')

%% FEEDBACK NORMALIZATION

