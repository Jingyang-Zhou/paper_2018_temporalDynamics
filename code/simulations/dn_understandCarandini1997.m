% dn_understandCarandini1997

%% The model:

% The driving current and the conductance vary over time

%% DEFINE STIMULUS AND A TIME COURSE

dt   = 0.001;
t    = dt : dt : 2;

% make a static 100ms stimulus (contrast = 1;)
stim = zeros(1, length(t));
stim(t>0.2 & t<=01.5) = 1;

%% PRE-DEFINED PARAMETERS

prm.k    = [0.1, 0.5, 0.8]; % effectiveness of the normalization pool, so should be between 0 and 1?
prm.g0   = 0.1; % a constant, baseline conductance?

prm.w    = 0.9; % determins the shape of the impulse repsonse function
prm.tau1 = 0.05; % summation length of the impulse response function
prm.C    = 2; % capacitance

% USEFUL FUNCTIONS --------------------------------------------------------
compute_g    = @(g0, k, R) g0./sqrt(1 - k .* R);
compute_g_CN = @(g0, k, L) g0 + k * L;

%% MAKE THE IMPULSE RESPONSE FUNCTION

irf1 = t.*exp(-t./prm.tau1); 
irf1 = irf1./sum(irf1);

irf2 = t .* exp(-t./(prm.tau1 * 1.5)); 
irf2 = irf2./sum(irf2);

irf = irf1 - prm.w * irf2;

figure (1), clf
subplot(2, 2, 1)
plot(t, irf, 'k-', 'linewidth', 2), box off, title('impulse response function')
xlabel('time (s)'), set(gca, 'fontsize', 14)

%% INITIATE THE MODEL PREDICTION

nTimeCourses = length(prm.k);

delta_V = 0;
V       = zeros(nTimeCourses, length(t)); % membrane potential
g       = zeros(nTimeCourses, length(t)); % conductance
R       = zeros(nTimeCourses, length(t)); % firing rate

%% COMPUTE MODEL PREDICTION

I = convCut(stim, irf, length(t));

for k = 1 : nTimeCourses
    for it = 1 : length(t) - 1
        g(k, it) = compute_g(prm.g0, prm.k(k), R(k, it));
        delta_V   = (I(it) - g(k, it) * V(k, it))/prm.C;
        V(k, it + 1) = V(k, it) + delta_V;
        R(k, it + 1) = max(V(k, it + 1),0)^2;
    end
end

R = R./max(R(:));
V = V./max(V(:));

%% VISUALIZE MODEL PREDICTIONS

figure (2), clf
subplot(3, 1, 1)
% plot stimulus
plot(t, stim, 'k-', 'linewidth', 2), box off, hold on
% plot membrane potential
plot(t, V, '-', 'linewidth', 2, 'color', 'r')
% plot spike rate
plot(t, R, '-', 'linewidth', 2, 'color', 0.2* ones(1, 3))
%legend('stimulus', 'membrane potential', '','spike rate')
set(gca, 'fontsize', 14)
title('stimulus, MP and spike rate')

subplot(3, 1, 2), set(gca, 'colororder', copper(nTimeCourses)), hold on
plot(V(:, t>0.2 & t<0.3)', 'linewidth', 3),
xlabel('time (ms)'), ylabel('Membrane potential'), legend('0.1', '0.5', '0.8', 'Location', 'best'), 
set(gca, 'fontsize', 14), title('Effect of conductance on MP'), xlim([0, 70])

subplot(3, 1, 3), 
scatter(R(3, :), g(3, :)), box off, 
xlabel('total firing rate'), ylabel('conductance')
set(gca, 'fontsize', 14)

%% COMPRESSIVE NONLINEARITY MODEL

nTimeCourses = length(prm.k);

delta_V = 0;
V       = zeros(nTimeCourses, length(t)); % membrane potential
g       = zeros(nTimeCourses, length(t)); % conductance
R       = zeros(nTimeCourses, length(t)); % firing rate

%% COMPUTE MODEL PREDICTION

% compute_g_CN

I = convCut(stim, irf, length(t));

for k = 1 : nTimeCourses
    for it = 1 : length(t) - 1
        g(k, it) = compute_g_CN(prm.g0, prm.k(k), 1);
        delta_V   = (I(it) - g(k, it) * V(k, it))/prm.C;
        V(k, it + 1) = V(k, it) + delta_V;
        R(k, it + 1) = max(V(k, it + 1),0)^2;
    end
end

R = R./max(R(:));
V = V./max(V(:));

figure (3), clf
subplot(3, 1, 1)
% plot stimulus
plot(t, stim, 'k-', 'linewidth', 2), box off, hold on
% plot membrane potential
plot(t, V, '-', 'linewidth', 2, 'color', 'r')
% plot spike rate
plot(t, R, '-', 'linewidth', 2, 'color', 0.2* ones(1, 3))
%legend('stimulus', 'membrane potential', '','spike rate')
set(gca, 'fontsize', 14)
title('stimulus, MP and spike rate')

subplot(3, 1, 2), set(gca, 'colororder', copper(nTimeCourses)), hold on
plot(V(:, t>0.2 & t<0.3)', 'linewidth', 3),
xlabel('time (ms)'), ylabel('Membrane potential'), legend('0.1', '0.5', '0.8', 'Location', 'best'), 
set(gca, 'fontsize', 14), title('Effect of conductance on MP')

subplot(3, 1, 3), 
scatter(R(3, :), g(3, :)), box off, 
xlabel('total firing rate'), ylabel('conductance')
set(gca, 'fontsize', 14)

%% COMPUTE THE MODEL PREDICTION TO STIMULUS OF DIFFERENT CONTRASTS
% 
contrast = linspace(0.01, 1, 10);
V       = zeros(length(contrast), length(t));
g       = zeros(length(contrast), length(t));
R       = zeros(length(contrast), length(t));

R_max = []; t_max = [];

for k = 1 : length(contrast)
    I = convCut(stim .* contrast(k), irf, length(t));
    delta_V = 0;
    for it = 1 : length(t) - 1
        g(k, it) = compute_g_CN(prm.g0, 0.7, R(k, it));
        delta_V   = (I(it) - g(k,it) * V(k,it))/prm.C;
        V(k,it + 1) = V(k,it) + delta_V;
        R(k, it + 1) = max(V(k, it + 1),0)^2;
    end
    R_max(k, :) = R(k, :)./max(R(k, :));
    t_max(k) = t(find(R(k, :) == max(R(k, :))));
end

%% VISUALIZE CONTRAST-DEPENDNET MODEL PREDICTIONS

figure (4), 
subplot(2, 2, 1), cla, set(gca, 'colorOrder', copper(length(contrast))), hold on
plot(t, R', 'linewidth', 2), title('Un-normalized response time courses')
xlim([0.2, 0.6]), legend('0.12', '0.25', '0.5', '1')
set(gca, 'fontsize', 14), xlim([0.2, 0.7])

subplot(2, 2, 2), cla, set(gca, 'colorOrder', copper(length(contrast))), hold on
plot(t, R_max', 'linewidth', 2)
xlim([0.2, 0.6]), title('Normalized response time courses')
set(gca, 'fontsize', 14), xlim([0.2, 0.7])

% plot response amplitude versus phase
subplot(2, 2, 3), cla
plot(contrast, sum(R'), 'ko-', 'linewidth', 2), box off
xlabel('contrast'), ylabel('maximal response'), set(gca, 'fontsize', 14)
title('contrast-dependent response level')

subplot(2, 2, 4),  cla
plot(contrast, fliplr(t_max), 'ko-', 'linewidth', 2), box off
xlabel('contrast'), ylabel('maximal response time (s)'), set(gca, 'fontsize', 14)
title('contrast-dependent max. response time'), 