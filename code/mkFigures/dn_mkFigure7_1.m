% I think this file aims to compare 3 models. Heeger 1992, Carandini 1997,
% and Mikaelian and Simoncelli 2001.

% Heeger model has feedback component, not of a RC circuit type
% The other two models

% I will compare the two RC-circuit models, then compared them back to the
% feedback implementation.

% Maybe I can test the steady state prediction to contrast for each model.

%% CHIRP STIMULUS?

t    = 0.001 : 0.001 : 10;
stim = zeros(size(t));

% first 5 seconds: a short and a long stimulus
stim(201 : 1100)  = 1;
stim(1801 : 1900) = 1;

% next, stimulus with brief gaps
stim(3501 : 4000) = 1;
stim(4301 : 4800) = 1;

% last, stimulus with different contrast
stim(6501 : 7000) = 0.05;
stim(7501 : 8000) = 0.2;
stim(8501 : 9000) = .5;

normMax = @(x) x ./max(x);

figure (1), subplot(8, 1, 1), cla, plot(t, stim, 'k-', 'linewidth', 2.5), box off, title('stimulus')
set(gca, 'fontsize', 14), ylabel('contrast')

%% INJECTED CURRENT

% Here, I assume that the current injected is a slightly low-passed
% stimulus squarewave. This is for taking account into the front end
% processing. May or may not be important.

% t    = 0.001 : 0.001 : 2;
% stim = zeros(size(t)); stim(1 : 600) = 1;stim(881 : 1000) = 1;

irf  = gammaPDF(t, 0.05, 2); % assuming very brief impulse
I    = convCut(stim, irf, length(t));

%% IMPLEMENT THE CARANDINI 97 MODEL

% Here the transform between V and R is assumed to be half-wave
% rectification
compute_g = @(g0, k, V) g0./sqrt(1 - k .* max(0, V));

g  = zeros(size(t));
V1 = 0.1 * ones(size(t));

k  = 1.2;
g0 = 1;
C  = 2; % capacitance

for it = 1 : length(t) - 1
    g(it)      = compute_g(g0, k, V1(it));
    dV         = (I(it) - g(it) * V1(it))/C;
    V1(it + 1) = V1(it) + dV;
end

figure (1), subplot(7, 1, 2), cla, plot(t, normMax(V1), 'r-', 'linewidth', 2.5), hold on, plot(t, normMax(I), 'k-'), axis tight, box off
title('Carandini and Heeger 1994'), set(gca, 'fontsize', 14)
% subplot(4, 2, 2), plot(t, g)

%% HEEGER 1992 MODEL

alpha = 0.01;
delta = 3;
R_max = 1; 
sigma = 0.1;

R_t = zeros(size(t)); B_t = zeros(size(t));

for it = 2 : length(t)
    B_t(it) = alpha * R_t(it) + (1- alpha) * B_t(max(it - delta, 1)); % this is where the exponential decay comes from
    R_t(it + delta) = I(it)./sigma * (R_max - B_t(it));
end

figure (1), subplot(7, 1, 3), cla, plot(t, normMax(R_t(1 : length(t))), 'r-', 'linewidth', 2.5), hold on, plot(t, normMax(I), 'k-'), axis tight, box off
title('Heeger 1992'), set(gca, 'fontsize', 14)

%% IMPLEMENT MIKAELIAN-SIMONCELLI MODEL

del_t = 0.03;
del_t = del_t ./(t(2) - t(1));

% weight of the response w here is the same as k in the previous model
w     = 0.01;
V2    = .1* (size(t));
g     = .1 * ones(size(t));
C     = 3;
% Assume half-wave rectification between V and R (in the original paper, a squaring after half-wave rectification is assumed)

for it = 1 : length(t) - 1
    g_pas     = g(max(it - del_t, 1));
    V_pas      = V2(max(it - del_t, 1));
    g(it)      = g_pas * w * max(V_pas, 0);
    g(it)      = sqrt(g(it));
    g(it)      = mean([g(max(1, it -15) : it)]);% low pass
    dV         = (I(it) - g(it) * V2(it))/C;
    V2(it + 1) = V2(it) + dV;
end
figure (1), subplot(7, 1, 4), cla, plot(t, normMax(V2), 'r-', 'linewidth', 2.5), hold on, plot(t, normMax(I), 'k-'), axis tight, box off
title('Mikaelian and Simoncelli 2001'), set(gca, 'fontsize', 14)
% subplot(4, 2, 4), plot(t, g)

%% IMPLETEMENT MANTE 2008

r = []; g= [];

Cc = 1; % conductance for the contrast gain control
r  = I;

g_scale = 1;
lp      = gammaPDF(t, 0.1, 2);
g       = g_scale * convCut(I, lp, length(t));

for it = 1 : length(t) -1
    del_r = 1./Cc.* (I(it) - g(it) * r(it));
    r(it+1) = r(it) + del_r;
end

figure (1), subplot(7, 1, 5), cla, plot(t, normMax(r), 'r-', 'linewidth', 2.5), hold on, plot(t, normMax(I), 'k-'),  box off, axis tight
title('Mante et al. 2008'), set(gca, 'fontsize', 14)
% subplot(4, 2, 6), plot(t, g),
%% 2TC MODEL

params = [];
params.b1 = 0.8;
params.b2 = 0.4;

pred = dn_2Chansmodel(params, stim, t, 0.001);

figure (1), subplot(7, 1, 6), cla, plot(t, normMax(pred), 'r-', 'linewidth', 2.5), hold on, plot(t, normMax(I), 'k-'),  box off, axis tight
title('Horiguchi and Wandell 2009'), set(gca, 'fontsize', 14)

%% DN MODEL

prm = [0.05, 0, 0.15, 2, 0.1, 0, 1];
dn = dn_DNmodel(prm, stim, t);

figure (1), subplot(7, 1, 7), cla, plot(t, normMax(dn), 'r-', 'linewidth', 2.5), hold on, plot(t, normMax(I), 'k-'),  box off, axis tight
title('DN model'),  set(gca, 'fontsize', 14)
%subplot(4, 2, 8), plot(t, g),
