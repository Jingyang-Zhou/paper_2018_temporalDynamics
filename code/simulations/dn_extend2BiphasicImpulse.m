% DN_simulations

% This script is meant to test the most sensible way to generalize the DN
% model to include biphase impulse response function. The current issue is
% that when we use biphasic impulse response function, we fix n to be 2 to
% avoid dealing with complex numbers. The more general way to do it is to
% put an absolute value somewhere, and this script is used to figure out
% where. 

%% PRE-DEFINE MODEL INPUTS AND PARAMETERS

% time course and stimulus
t    = 0.001 : 0.001 : 1.2;
stim = zeros(size(t));
stim(201 : 700) = 1;

% DN model parameters
prm = [];
prm.tau1  = 0.05;
prm.w     = 0.9;
prm.tau2  = 0.1;
prm.n     = 1.5;
prm.sigma = 0.1;

% useful function
normx = @(x) x./norm(x, 1);

%% DERIVED PARAMETERS

% impulse response function
irf = gammaPDF(t, prm.tau1, 2) - prm.w * gammaPDF(t, prm.tau1 * 1.5, 2);

% delayed filter
del = normx(exp(-t/prm.tau2));

% visualize
figure (1), clf
subplot(4, 2, 1), plot(t, irf, 'k-', 'linewidth', 3), axis tight, box off, xlabel('t (s)'), title('IRF')

%% COMPUTE DIFFERENT STAGES OF MODEL PREDICTION

rsp1 = [];

component1 = convCut(irf, stim, length(irf));
component2 = convCut(component1, del, length(irf));

% NUMERATOR
rsp1.num = component1.^prm.n;
% DENOMINATOR
rsp1.den = prm.sigma^prm.n + component2.^prm.n;
% NORMALIZED RESPONSE
rsp1.norm = abs(rsp1.num ./rsp1.den);

subplot(4, 2, 3), plot(t, rsp1.num)
subplot(4, 2, 5), plot(t, rsp1.den)
subplot(4, 2, 7), plot(t, rsp1.norm)

rsp2 = [];

% NUMERATOR
rsp2.num = abs(component1.^prm.n);
% DENOMINATOR
rsp2.den = prm.sigma^prm.n + abs(component2.^prm.n);
% NORMALIZED RESPONSE
rsp2.norm = rsp2.num ./rsp2.den;

subplot(4, 2, 4), plot(t, rsp2.num), title('numerator')
subplot(4, 2, 6), plot(t, rsp2.den), title('denominator')
subplot(4, 2, 8), plot(t, rsp2.norm), title('final response')



