% dn_chooseElectrodes

%% SAVE FIGURE AND DATA KNOB

saveFigure = 1;
saveData   = 0;

%% PREDEFINED VARIABLES

prm.ecog.roiNm = {'V1', 'V2', 'V3', 'lateral', 'ventral', 'dorsal'};
prm.ecog.selThresh = 0.5;
prm.ecog.tBase     = 1 : 200;

%% USEFUL FUNCTIONS

normBase = @(x, range) x - mean(x(range));
normMax  = @(x) x./max(x);

%% LOAD DATA

dataLoc = fullfile(dn_ECoG_RootPath, 'data');
fName   = 'dn_preprocessedData.mat';

a    = load(fullfile(dataLoc, fName));
dt   = a.dt;
data = dt.ecog;

prmFileNm = 'dn_params.mat';
b    = load(fullfile(dataLoc, prmFileNm));
prm  = b.prm;

%% DERIVED VARIABLES

label  = data.labels;
elec   = data.chans;
bb     = data.bb;
stimNm = data.stimNm;

t_lth  = size(bb{1}, 2);
nElec  = size(bb{1}, 3);

prm.ecog.srate = 1000;
prm.ecog.t     = 1/prm.ecog.srate : 1/prm.ecog.srate : t_lth/prm.ecog.srate;
prm.ecog.stim  = zeros(1, t_lth);
prm.ecog.stim(201 : 700) = 1;

roiNms = {'V1', 'V2', 'V3', 'V3A', 'hV4', 'VO', 'LO', 'TO', 'IPS'};
nrois  = length(roiNms);
roiIdx = {1, 2, 3, 7, [5, 6], [4, 8]};

combined_nrois = length(roiIdx);

%% PREP TO CHANGE TIME COURSE BROADBAND SIGNAL TO PERCENTAGE CHANGE

% Here, we only look at response to noise stimuli. For analysis comparing
% broadband repsonse to grating versus broadband response to noise, see
% (1) dn_comparebb_GratingNoise.m
% (2) dn_spectrogram_spectra.m

mbb{1} = squeeze(mean(bb{1}(stimNm < 4, :, :))); % Two different frequency ranges
mbb{2} = squeeze(mean(bb{2}(stimNm < 4, :, :)));

% VISUALIZE "mbb" ---------------------------------------------------------
% figure
% for iElec = 1 : nElec
%    subplot(8, 10, iElec)
%    plot(mbb{1}(:, iElec), 'r-', 'linewidth', 2), hold on
%    plot(mbb{2}(:, iElec), 'k-', 'linewidth', 2), xlim([0, 1204]), box off
% end
% WHAT I LEARNED FROM THIS PLOT -------------------------------------------
% (1) The bb time courses using two ranges of frequency bands have
%     different baseline and maximal response, although they look similar in
%     shape.
% (2) I want to compare the results of bb extracted using two frequency
%     ranges but with the same electrodes. Using different electrodes for
%     different bb extraction condition introduces unnecessary
%     complications in the analysis.
% (3) From now on, in this script, I will only use the first way of
%     extracting broadband, which is with frequency range [70, 210].

%% CHANGE TIME COURSES TO PERCENTAGE SIGNAL AND SELECT ELECTRODE

% THE QUESTION IS AS THE FOLLOWING:
% If there are two electrodes, the first one has a baseline level around 2,
% and the stimulus-triggered peak response level around 4. The second
% electrode has a baseline level 20, and peak response level 22. Should we
% deduce baseline levels from both electrodes and treat them as having the
% same response, or should be compare only the percentage change level so
% that the first electrode is treated as having a better signal than the
% second?

% I think the second way should be preferred. Here is the reasoning:
% imagine different baseline level represents different resting states at
% differet parts of the brain, some areas have higher resting state,
% therefore higher baseline, and possibly higher response variance at the
% baseline as well. Therefore the same absolute value of signal increase
% given high baseline level could be because of noise, while given low
% baseline level could be due to a big signal chance. This is perhaps not
% the best argument, need to improve later.

select_idx = [];

for iElec = 1 : nElec
    % CONVERT THE SIGNAL TO PERCENTAGE CHANGE FROM THE BASELINE
    baseline = mbb{1}(prm.ecog.tBase, iElec);
    m_baseline(iElec) = mean(baseline);
    s_baseline(iElec) = std(baseline);
    
    mbb{1}(:, iElec) = mbb{1}(:, iElec)/m_baseline(iElec);
    % SUBTRACT THE BASELINE
    mbb{1}(:, iElec) = mbb{1}(:, iElec) - 1;
    
    % ALTERNATIVELY (OBSOLETE), USE ABSOLUTE SIGNAL RATHER THAN %CHANGE
    %  mbb{1}(:, iElec) = mbb{1}(:, iElec) - baseline;
end

%% SELECT ELECTRODE BASED ON THE MAX OF THE STIMULUS-TRIGGERED PERCENTAGE CHANGE SIGNAL

% FIRST SELECTION CRITERION:
select_idx1 = max(mbb{1}) > prm.ecog.selThresh;

% SECOND SELECTION CRITERION:
stim_on = 200 : 700;

for k = 1 : length(elec)
    if mean(mbb{1}(stim_on, k)) > 0,
        select_idx2(k) = 1;
    else
        select_idx2(k) = 0;
    end
end

select_idx = select_idx1 + select_idx2;
select_idx = (select_idx == 2);

%% JUSTIFYING WHY CHOOSE PERCENTAGE CHANGE FROM THE BASELINE INSTEAD OF ABSOLUTE CHANGE

% If baseline std. is independent of baseline mean, then we should use
% absolute change, if baseline std. increases with baseline mean, then we
% should use percentage change. This by no means be the only reason why
% should or shouldn't we use the percentage change signal.

figure (100), clf
scatter(m_baseline, s_baseline, 'ko', 'markerfacecolor', 'r')
xlabel('baseline mean'), ylabel('baseline std.'),
title('Justify why we use percentage change signal'), set(gca, 'fontsize', 14)

%% VISUALIZE THE BASELINE NORMALIZED "mbb"

fg2 = figure (2); clf, fg2.Position = [1, 2000, 2000, 2000];
for iElec = 1 : nElec
    subplot(8, 10, iElec),
    plot(mbb{1}(:, iElec), 'k-', 'linewidth', 2), xlim([1, 1204]), hold on,
    plot([1, 1204], [prm.ecog.selThresh, prm.ecog.selThresh], 'r-')
    ylim([-2, 10]), box off, set(gca, 'xtick', ''), title(data.chans(iElec))
    if select_idx(iElec) == 1, set(gca, 'color', 'y'), end, set(gca, 'fontsize', 12)
end

%% MODIFY ELECTRODE NAMES

% Hand-modify the electrodes the names of which do not match the cortical
% locations. (VO and hV4 labels)
chans_mod = [103, 104, 99, 100, 94, 95, 57, 63,  62,  59, 105, 117];

for k = 1: length(chans_mod)
    modidx(k) = find(elec == chans_mod(k));
end

labs_mod = {'VO', 'VO', 'hV4', 'hV4', 'hV4', 'hV4', 'V3', 'LO', 'LO', 'LO',  '', ''};
label(modidx) = labs_mod;

%% COMPUTE SELECTED LABELS AND ELECTRODES, TIME COURSES

select_Labels   = label(select_idx);
select_bbts_tmp = mbb{1}(:, select_idx);
elec_idx        = elec(select_idx);

%% LABEL SELECTED ELECTRODES WITH ROI NAMES

nm_idx = {};

for iroi = 1 : nrois
    nm_idx{iroi} = find(not(cellfun('isempty', strfind(select_Labels,roiNms{iroi}))))
    if iroi == 3,
        v3aIdx       = find(not(cellfun('isempty', strfind(select_Labels,roiNms{4}))));
        intersectIdx = find(not(cellfun('isempty', strfind(select_Labels,'V3-p/V3A'))));
        nm_idx{iroi} = [setdiff(nm_idx{iroi}, v3aIdx), intersectIdx]
    end;
end

% COMBINE INDICES BECAUSE WE DONT' HAVE ENOUGH ELECTRODES IN THE MORE ANTERIOR REGIONS
% roiIdx = {1, 2, 3, 7, [5, 6], [4, 8]};
cb_nmIdx = {};
cb_nmIdx(1:3) = nm_idx(1:3);
cb_nmIdx(4)   = nm_idx(7);
cb_nmIdx{5}   = [nm_idx{5}, nm_idx{6}];
cb_nmIdx{6}   = [nm_idx{4}, nm_idx{8}];

% DISPLAY ROI NAMES FOR THE SAKE OF CHECKING ------------------------------
display('[dn_chooseElectrodes]: Electrode labels in each ROI:')
for iroi = 1 : combined_nrois
    disp(prm.ecog.roiNm{iroi})
    disp(select_Labels(cb_nmIdx{iroi}))
end
%% COMBINE ROI BROADBAND TIME COURSES AND ELECTRODE INDICES

final_bbts = {};

for iroi = 1 : length(cb_nmIdx)
    final_bbts{iroi} = select_bbts_tmp(:, cb_nmIdx{iroi});
    final_elec{iroi} = elec_idx(cb_nmIdx{iroi});
end

%% VISUALIZE FINAL SELECTED BROADBAND TIME SERIES

% NORMALIZE EACH TIME COURSE TO THE MAX -----------------------------------
for iroi = 1 : length(cb_nmIdx)
    for k = 1 : size(final_bbts{iroi}, 2)
        nm_final_bbts{iroi}(:, k) = normMax(final_bbts{iroi}(:, k));
    end
end

figure (3), clf
for iroi = 1 : length(cb_nmIdx)
    subplot(2, 3, iroi)
    to_plot = normMax(mean(nm_final_bbts{iroi}, 2));
    plot(to_plot, 'k-', 'linewidth', 2), xlim([0, 1204]), hold on
    if iroi == 1,
        max_v1 = find(to_plot ==max(to_plot)),
        off_v1 = to_plot(700);
    end
    plot([max_v1, max_v1], [0, max(to_plot)], 'r-'),
    plot([0, 1204], [[off_v1], off_v1], 'r-')
    plot([700, 700], [0, 1], 'r-'), ylim([-0.1, 1])
    set(gca, 'xaxislocation', 'origin'), box off
end

%% SAVE FIGURES

if saveFigure
    figLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
    
    % save figure 2: the electrode selection figure
    fg2Nm  = 'pre_electrodeSelection';
    printnice(2, 0, figLoc, fg2Nm);
    % save the figure that shows the signal to noise ratio in the baseline
    fg100Nm = 'pre_sig2noise_baseline';
    printnice(100, 0, figLoc, fg100Nm);
end

%% SAVE DATA
% I want to save the data (pre-normalized), corresponding electrode labels,
% and the prm
if saveData
    % SAVE DATA -----------------------------------------------------------
    dt.ecog = data;
    dt.ecog.bbts_roi = final_bbts;
    dt.ecog.elec_roi = final_elec;
    save(fullfile(dataLoc, fName), 'dt');
    % SAVE PARAMS ---------------------------------------------------------
    prmFileNm = 'dn_params.mat';
    save(fullfile(dataLoc, prmFileNm), 'prm')
end
