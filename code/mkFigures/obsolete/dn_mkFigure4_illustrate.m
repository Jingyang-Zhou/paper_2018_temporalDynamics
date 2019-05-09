% dn_mkFigure4_illustration

%% MAKE ILLUSTRATION FOR INCREASING TIME TO PEAK

dt = 0.001;
t  = dt : dt : 4;
stim = zeros(size(t));
stim (200 : end) = 1;

% model parameteres 
tau1 = [0.2, 0.05];
sigma = [0.04, 0.1];

figure (1), clf
for k = 1 : 2
    dnPrm = [tau1(k), 0, 0.1, 2, sigma(k), 0, 1];
    dnPrd = normMax(dn_DNmodel(dnPrm, stim, t));
    subplot(2, 1, k)
    plot(t, dnPrd, 'r-', 'linewidth', 3), hold on, box off, 
    xlim([0.18, 0.4]), axis off
end

saveLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
printnice(1, 0, saveLoc, 'nfg_figure4A1')

%% MAKE ILLUSTRATION FOR ASYMPTOTIC responses

% model parameteres 
sigma = [0.3, 0.005];
tau1  = [0.1, 0.2];

figure (2), clf
for k = 1 : 2
    dnPrm = [tau1(k), 0, 0.1, 2, sigma(k), 0, 1];
    dnPrd = normMax(dn_DNmodel(dnPrm, stim, t));
    subplot(2, 1, k)
    plot(t, dnPrd, 'r-', 'linewidth', 3), hold on, box off, 
    xlim([0.18, 0.8]), axis off
end

saveLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
printnice(2, 0, saveLoc, 'nfg_figure4A2')

%% MAKE ILLUSTRATION FOR IMPULSE RESPONSE FUNCTION SHAPE

figure (3), clf
tau1 = 0.2;
irf(1, :) = normMax(gammaPDF(t, tau1, 2));
irf(2, :) = normMax(gammaPDF(t, tau1, 2) - gammaPDF(t, 0.3, 2));

for k = 1 : 2
   subplot(2, 1, k)
   plot(t, irf(k, :), 'g-', 'linewidth', 3), box off, xlim([0, 2]), axis off
end
saveLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
printnice(3, 0, saveLoc, 'nfg_figure4B')
