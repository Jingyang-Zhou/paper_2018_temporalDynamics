% DN_mkNewFigure2

saveFigure = 0;

%% MAKE IMPULSE RESPONSE FUNCTION

t = 0.001 : 0.001 : 1;

irf1   = gammaPDF(t, 0.05, 2);
irf1_1 = gammaPDF(t, 0.1, 2);
irf2   = irf1 - 0.8 * irf1_1;

figure (1), clf
subplot(3, 2, 1), plot(t, irf1, 'k-', 'linewidth', 3)
subplot(3, 2, 2), plot(t, irf2, 'k-', 'linewidth', 3),
for k = 1 : 2
    subplot(3, 2, k)
    axis tight, box off, axis off, axis square
end

%% COMPUTE LINEAR REPSONSE WITH BOTH IMPULSE REPSONSE FUNCTION TYPE

stim = zeros(size(t));

stim (101 : 500) = 1;

lin = [];
lin(1, :) = normMax(convCut(irf1, stim, length(irf1)));
lin(2, :) = normMax(convCut(irf2, stim, length(irf1)));

figure (1),
for k = 1 : 2
    subplot(3, 2, k + 2)
    plot(t, lin(k, :), 'k-', 'linewidth', 3),
    axis tight, ylim([-0.6, 1]), box off, axis square, axis off
end

%% COMPUTE RECTIFIED RESPONSE

figure (1),
for k = 1 : 2
    subplot(3, 2, k + 4)
    plot(t, max(lin(k, :), 0), 'k-', 'linewidth', 3),
    axis tight, ylim([-0.6, 1]), box off, axis square, axis off
end

%% FIGURE 2

param = [0.05, 0, 0.1, 2, 0.1, 0, 1];

prd(1, :) = normMax(dn_DNmodel(param, stim, t));

newParam = param; newParam(1) = 0.15;
prd(2, :) = normMax(dn_DNmodel(newParam, stim, t));

newParam = param; newParam(4) = 3;
prd(3, :) = normMax(dn_DNmodel(newParam, stim, t));

newParam = param; newParam(5) = 0.01;
prd(4, :) = normMax(dn_DNmodel(newParam, stim, t));

newParam = param; newParam(3) = 0.5;
prd(5, :) = normMax(dn_DNmodel(newParam, stim, t));

figure (2), clf
for k = 1 : 4
    subplot(2, 4, k),
    plot(t, stim, 'k-'), hold on, plot(t, prd(1, :), 'r-', 'linewidth', 3), axis square, box off, axis off, ylim([0, 2])
end
subplot(2, 4, 5)
plot(t, stim, 'k-'), hold on, plot(t, prd(2, :), 'r-', 'linewidth', 3), axis square, box off, axis off, ylim([0, 2])

subplot(2, 4, 6)
plot(t, stim, 'k-'), hold on, plot(t, prd(3, :), 'r-', 'linewidth', 3), axis square, box off, axis off, ylim([0, 2])

subplot(2, 4, 7)
plot(t, stim, 'k-'), hold on, plot(t, prd(4, :), 'r-', 'linewidth', 3), axis square, box off, axis off, ylim([0, 2])

subplot(2, 4, 8)
plot(t, stim, 'k-'), hold on, plot(t, prd(5, :), 'r-', 'linewidth', 3), axis square, box off, axis off, ylim([0, 2])

%%

if saveFigure
    saveLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
    fig1Nm  = 'nfg_figure2A';
    printnice(1, 0, saveLoc, fig1Nm);
    
    fig2Nm  = 'nfg_figure2B';
    printnice(2, 0, saveLoc, fig2Nm);
end