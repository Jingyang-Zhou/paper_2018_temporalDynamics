% dn_generateNewFigure1

saveFigure = 0;

%% PRE-DEFINED MODEL PARAMETERS

params = [0.1, 0, 0.1, 2, 0.1, 0, 1];
t      = 0.001 : 0.001 : 1.5;

normMax = @(x) x./max(x);

%% MAKE STIMULUS

stim = [];

% one 100ms pulse ---------------------------------------------------------
stim.onepulse_100 = zeros(1, length(t));
stim.onepulse_100(t > 0.1 & t <= 0.2) = 1;

% one 200ms pulse ---------------------------------------------------------
stim.onepulse_200 = zeros(1, length(t));
stim.onepulse_200(t > 0.1 & t <= 0.3) = 1;

% one sustained pulse -----------------------------------------------------
stim.onepulse_sustained = zeros(1, length(t));
stim.onepulse_sustained(t > 0.1 & t <= 0.8) = 1;

% two 100ms pulses, short ISI ---------------------------------------------
stim.twopulses_100 = zeros(1, length(t));
stim.twopulses_100(t > 0.1 & t <= 0.2) = 1;
stim.twopulses_100(t > 0.3 & t <= 0.4) = 1;

% two 100ms pulses, long ISI ---------------------------------------------
stim.twopulses_800 = zeros(1, length(t));
stim.twopulses_800(t > 0.1 & t <= 0.2) = 1;
stim.twopulses_800(t > 1 & t <= 1.1) = 1;

% one sustained pulse, high contrast --------------------------------------
stim.highContrast = zeros(1, length(t));
stim.highContrast(t > 0.1 & t <=0.8) = 1;

% one sustained pulse, low contrast ---------------------------------------
contrast = 0.3;
stim.lowContrast = stim.highContrast .* contrast ;

%% VISUALIZE STIMULI ------------------------------------------------------

figure (1), clf
subplot(4, 2, 1), plot(t, stim.onepulse_100, 'k-'),
subplot(4, 2, 2), plot(t, stim.onepulse_sustained, 'k-'),

subplot(4, 2, 3), plot(t, stim.onepulse_100, 'k-'),
subplot(4, 2, 4), plot(t, stim.onepulse_200, 'k-')

subplot(4, 2, 5), plot(t, stim.twopulses_800, 'k-')
subplot(4, 2, 6), plot(t, stim.twopulses_100, 'k-')

subplot(4, 2, 7), plot(t, stim.highContrast, 'k-')
subplot(4, 2, 8), plot(t, stim.lowContrast, 'k-')

for k = 1 : 8
    subplot(4, 2, k), hold on, ylim([0, 1]), box off, 
end

%% MAKE LINEAR PREDICTIONS

irf = gammaPDF(t, 0.05, 2);

lin = [];

% FIRST ROW ---------------------------------------------------------------
lin.onepulse_100       = convCut(irf, stim.onepulse_100, length(t));
lin.onepulse_sustained = convCut(irf, stim.onepulse_sustained, length(t));

maxSustainedRsp        = max(lin.onepulse_sustained);
lin.onepulse_100       = lin.onepulse_100./maxSustainedRsp  ;
lin.onepulse_sustained = normMax(lin.onepulse_sustained);

% SECOND ROW --------------------------------------------------------------
lin.onepulse_200      = convCut(irf, stim.onepulse_200, length(t));
lin.onepulse_200      = lin.onepulse_200 ./ maxSustainedRsp;

% THIRD ROW ---------------------------------------------------------------
lin.twopulses_800     = convCut(irf, stim.twopulses_800, length(t));
lin.twopulses_800     = lin.twopulses_800./maxSustainedRsp;

lin.twopulses_100     = convCut(irf, stim.twopulses_100, length(t));
lin.twopulses_100     = lin.twopulses_100./maxSustainedRsp;

%% VISUALIZE LINEAR PREDICTIONS

figure (1)
subplot(4, 2, 1), plot(t, lin.onepulse_100, 'r-', 'linewidth', 2)
subplot(4, 2, 2), plot(t, lin.onepulse_sustained, 'k:', 'linewidth', 2)

subplot(4, 2, 3), plot(t, lin.onepulse_100, 'r-', 'linewidth', 2)
subplot(4, 2, 4), plot(t, lin.onepulse_200, 'k:', 'linewidth', 2)

subplot(4, 2, 5), plot(t, lin.twopulses_800, 'r-', 'linewidth', 2)
subplot(4, 2, 6), plot(t, lin.twopulses_100, 'k:', 'linewidth', 2)

%% COMPUTE NON-LINEAR RESPONSE

% FIRST ROW ---------------------------------------------------------------
dn.onepulse_sustained = dn_DNmodel(params, stim.onepulse_sustained, t);
maxDNRsp = max(dn.onepulse_sustained);
dn.onepulse_sustained = normMax(dn.onepulse_sustained);

% SECOND ROW --------------------------------------------------------------
dn.onepulse_200 = dn_DNmodel(params, stim.onepulse_200, t);
dn.onepulse_200 = dn.onepulse_200./maxDNRsp;

% THIRD ROW ---------------------------------------------------------------
dn.twopulses_100 = dn_DNmodel(params, stim.twopulses_100, t);
dn.twopulses_100 = dn.twopulses_100./maxDNRsp;

% FOURTH ROW --------------------------------------------------------------
dn.highContrast   = dn_DNmodel(params, stim.highContrast, t);
dn.highContrast   = dn.highContrast./maxDNRsp;

dn.lowContrast   = dn_DNmodel(params, stim.lowContrast, t);
dn.lowContrast   = dn.lowContrast./maxDNRsp;

%% VISUALIZE NON-LINEAR RESPONSE

figure (1), 
subplot(4, 2, 2), plot(t, dn.onepulse_sustained, 'r-', 'linewidth', 2), title('Response declines for sustained stimulus')
subplot(4, 2, 4), plot(t, dn.onepulse_200, 'r-', 'linewidth', 2), title('subaddtive temporal sum')
subplot(4, 2, 6), plot(t, dn.twopulses_100, 'r-', 'linewidth', 2), title('adaptation')

% linear response
subplot(4, 2, 7), plot(t, dn.highContrast, 'r-', 'linewidth', 2)
subplot(4, 2, 8), plot(t, dn.highContrast .* contrast, 'k:', 'linewidth', 2)
subplot(4, 2, 8), plot(t, dn.lowContrast, 'r-', 'linewidth', 2), title('Contrast-dependent temporal dynamics')

%% SAVE FIGURE

if saveFigure,
    figLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
    fg1Nm = 'nfg_figure1A';
    printnice(1, 0, figLoc, fg1Nm);
    
    fg2Nm = 'nfg_figure1B';
    printnice(2, 0, figLoc, fg2Nm);
end

