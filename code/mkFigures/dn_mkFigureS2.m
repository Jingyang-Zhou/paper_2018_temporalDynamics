% dn_mkFigure_S4

saveFigure = 0;

%% LOAD DATA

dataFile = fullfile(dn_ECoG_RootPath, 'data', 'dn_preprocessedData.mat');
a = load(dataFile);

ecogdt = a.dt.ecog.bbts_roi;
elecNm = a.dt.ecog.elec_roi; 

%% PRE-DEFINED PARAMETERS

nrois = 4;
t     = [1 : 1204]./1000;

%% VISUALIZE RESPONSE FROM INDIVIDUAL ELECTRODE

fg1 = figure (1); clf
for iroi = 1 : nrois
    for k = 1 : size(ecogdt{iroi}, 2)
       subplot(4, 15, (iroi - 1) * 15 + k)
       patch([0.2, 0.7, 0.7, 0.2], [0, 0, 1, 1], 0.9  * ones(1, 3)), hold on
       plot(t, normMax(ecogdt{iroi}(:, k)), 'k-', 'linewidth', 3), axis tight, ylim([-0.4, 1]), box off
       set(gca, 'xaxislocation', 'origin'), set(gca, 'xtick', [0.2, 0.7], 'xticklabel', '', 'ytick', [0, 1], 'yticklabel', '')
       title(elecNm{iroi}(k)), set(gca, 'fontsize', 16)
    end
end

fg1.Position = [400, 2000, 2000, 400];


%% WITHOUT NORMALIZING TO THE MAX FOR EACH ELECTRODE

% fg2 = figure (2); clf
% for iroi = 1 : nrois
%     for k = 1 : size(ecogdt{iroi}, 2)
%        subplot(4, 15, (iroi - 1) * 15 + k)
%        patch([0.2, 0.7, 0.7, 0.2], [0, 0, 14, 14], 0.9  * ones(1, 3)), hold on
%        plot(t, ecogdt{iroi}(:, k), 'k-', 'linewidth', 3), axis tight, ylim([-0.4, 14]), box off
%        set(gca, 'xaxislocation', 'origin'), set(gca, 'xtick', [0.2, 0.7], 'xticklabel', '', 'ytick', [0, 14])
%        title(elecNm{iroi}(k)), set(gca, 'fontsize', 16)
%     end
% end
% 
% fg2.Position = [400, 2000, 2000, 400];


%% SAVE FIGURE

if saveFigure, 
    fg1Nm = 'fg_individualElectrode';
    fgLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
    printnice(1, 0, fgLoc, fg1Nm);
    
    fg2Nm = 'fg_individualElectrode_noNormalization';
    fgLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
    printnice(2, 0, fgLoc, fg2Nm);
end
    