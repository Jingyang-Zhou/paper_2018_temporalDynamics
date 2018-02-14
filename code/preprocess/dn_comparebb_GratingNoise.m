% dn_compareGratingNoise
%
% Here, we compare the broadband time course averaged across different
% types of stimuli. The prediction is that grating-triggered broadband time courses
% differ from the noise-triggered broadband time courses.

%% SAVE FIGURE KNOB

saveFigure = 1;

%% LOAD PRE-PROCESSED BROADBAND DATA

fName = 'dn_preprocessedData.mat';
fLoc  = fullfile(dn_ECoG_RootPath, 'data');
a     = load(fullfile(fLoc, fName));
data  = a.dt.ecog;

%% DERIVED VARIABLES

bb     = data.bb;
stimNm = data.stimNm;
nElec  = length(data.labels);

%% CATEGORIZE STIMULI

% We can tegorize stimuli in two ways:
% (1) 7 types, each stimulus as its own type
% (2) 2 types, noise and grating

bb_type  = {};

% THE FIRST TYPE OF STIMULUS GROUPING -------------------------------------
for iband = 1 : 2
    for iType = 1 : 7
        bb_type{iband, 1}(iType, :, :) = squeeze(mean(bb{iband}(stimNm == iType, :, :)));
    end
end
% THE SECOND TYPE OF STIMULUS GROUPING ------------------------------------
for iband = 1 : 2
    bb_type{iband, 2}(1, :, :) = squeeze(mean(bb{iband}(stimNm < 4, :, :)));
    bb_type{iband, 2}(2, :, :) = squeeze(mean(bb{iband}(stimNm > 3, :, :)));
end

%% VISUALIZE BROADBAND SIGNAL AGERAGED OVER DIFFERENT STIMULUS TYPES

% VISUALIZE BROADBAND SIGNAL FOR EACH INDIVIDUAL STIMULUS -----------------
fg1 = figure (1); clf, fg1.Position = [1, 2000, 2000, 2000];
for iElec = 1 : nElec
    subplot(8, 10, iElec), set(gca, 'colororder', copper(7)), hold on
    plot(squeeze(bb_type{1, 1}(:, :, iElec))', 'linewidth', 1.5), xlim([0, 1204])
    set(gca, 'xtick', [200, 700], 'xticklabel', ''), title(data.chans(iElec)), set(gca, 'fontsize', 12)
end
% VISUALIZE BROADBAND SIGNAL FOR 2 STIMULUS CLASSES -----------------------
fg2 = figure (2); clf, fg2.Position = [1, 2000, 2000, 2000];
for iElec = 1 : nElec
    subplot(8, 10, iElec), set(gca, 'colororder', copper(2)), hold on
    plot(squeeze(bb_type{1, 2}(:, :, iElec))', 'linewidth', 2), xlim([0, 1204])
    set(gca, 'xtick', [200, 700], 'xticklabel', ''), title(data.chans(iElec)), set(gca, 'fontsize', 12)
    if iElec == 1, legend('noise', 'grating'), end
end

%% SAVE FIGURES

if saveFigure
   fgNm1 = 'pre_bb_indiElec_stimClass_all';
   fgNm2 = 'pre_bb_indiElec_stimClass_noiseGrating';
   saveLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
   
   printnice(1, 0, saveLoc, fgNm1);
   printnice(2, 0, saveLoc, fgNm2);
end