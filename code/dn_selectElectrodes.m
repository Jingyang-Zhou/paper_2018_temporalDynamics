% Select electrodes for ECog experiment 1

%% DESCRIPTION

% This script does two things: select the electrodes within each ROI
% based on some criterion (maximum response of the baseline normalized
% broadband time course), as well as combining electrodes from 8 to 6
% different ROIs.

%% variables to be adjusted

thresh  = 3; % selection criterion, may want to update this to the data file

% equalize the baseline
t_preStim = 1 : 200;

%% load data

dataLoc = fullfile(dn_ECoG_RootPath, 'data');
fName   = 'dn_data.mat';
a       = load(fullfile(dataLoc, fName));
raw     = a.raw; 
% raw broadband data

%% Extract broadband parameters

% broadband data
goodElecs = raw.goodChannels; % 71 out of 118 electrodes
goodLabs  = raw.goodLabels; % labels for the 71 good electrodes

roiNm   = a.param.indiRoiNm; % roi names
nBoots  = 100; % number of bootstraps
nrois   = length(roiNm); % total number of ROIs

% averaged the data across trials
mdata = squeeze(mean(raw.bbmatrix, 2));
       
%% Select electrodes based on a threshold criterion

% Selection criterion: 
% First we equalize the baseline in each electrode, then we select the
% electrodes the maximum reposnse of which is above certain threshold

for k = 1 : size(mdata, 2)
    mdata_base0(:, k) = mdata(:, k) - mean(mdata(t_preStim, k));
end

thresh_idx  = max(mdata_base0) > thresh;

% extract corresponding data labels
selectLabs = goodLabs(thresh_idx); % selected labels
selectElec = goodElecs(thresh_idx); % selected electrode numbers
selectTS   = mdata_base0(:, thresh_idx); % bb time series of the selected electrodes

%% Label selected electrodes with ROI names

idx = []; 

for iroi = 1 : nrois
    idx{iroi} = find(not(cellfun('isempty', strfind(selectLabs,roiNm{iroi}))));
    if iroi ==3,
        v3aIdx       = find(not(cellfun('isempty', strfind(selectLabs,roiNm{4}))));
        intersectIdx = find(not(cellfun('isempty', strfind(selectLabs,'V3-p/V3A'))));
        idx{iroi}    = [setdiff(idx{iroi}, v3aIdx), intersectIdx];
    end
    % roi Labels
    indi_roiLabels{iroi} = selectLabs(idx{iroi});
    % time courses
    indi_roiTS{iroi} = selectTS(:, idx{iroi});
    % refer to the original 112 electrode index
end

%% combine ROIs

newRoiNm  = {'V1', 'V2', 'V3', 'lateral', 'ventral', 'dorsal'};
newRoiIdx = {1, 2, 3, 7, [5, 6], [4, 8]};

for k = 1 : length(newRoiIdx)
    if length(newRoiIdx{k}) == 1
        roiLabels{k} = indi_roiLabels{newRoiIdx{k}};
        roiTS{k}     = indi_roiTS{newRoiIdx{k}};
    elseif length(newRoiIdx{k}) == 2
        roiLabels{k} = [indi_roiLabels{newRoiIdx{k}(1)}, indi_roiLabels{newRoiIdx{k}(2)}];
        roiTS{k}     = [indi_roiTS{newRoiIdx{k}(1)}, indi_roiTS{newRoiIdx{k}(2)}];
    end
end

%% Bootstrap time series

bsData = [];

for k = 1 : length(newRoiNm)
    % generate indexes for bootstrapping
    nRoiElec = size(roiTS{k}, 2);
    bsidx    = randi([1, nRoiElec], nRoiElec, nBoots);
    for iBoot = 1 : nBoots
        ts       = mean(roiTS{k}(:, bsidx(:, iBoot)), 2);
        baseline = mean(ts(t_preStim));
        bsData(k, iBoot, :) = ts - baseline;
    end
end

%% visualize data

% to compare time series in later visual areas to V1
v1ts = mean(bsData(1, :, :), 2);
peak = find(v1ts == max(v1ts));
amp  = v1ts(700);

figure (1), clf

for k = 1 : length(newRoiNm)
    subplot(2, 3, k),
    m = squeeze(mean(bsData(k, :, :), 2));
    s = squeeze(std(bsData(k, :, :), [], 2));
    plot(m, 'k-'), hold on
    plot(m + s, 'color', .7*[1, 1, 1]);
    plot(m - s, 'color', .7*[1, 1, 1]);
    
    plot([peak, peak], [0, 15], 'r:'), plot([0, 1204], [0, 0], 'k:')
    plot([1, 1200], [amp, amp], 'r:'), plot([700, 700], [0, amp], 'r:')
    
    ax = gca;
    ax.XAxisLocation = 'origin'; box off
    set(gca, 'YLim', [-2 18]), xlim([0, 1204])
    title(newRoiNm{k})
end


%% re-save the data with only the relevant stimulus types

% param = [];
% param.thresh = thresh;
% param.roiNm  = newRoiNm;
% param.elecs      = selectElec;
% param.elecLabels = selectLabs;
% param.elecTS     = selectTS;
% param.roiLabels  = roiLabels;
% param.roiTS      = roiTS;
% 
% save(fullfile(dataLoc, fName), 'bsData', 'param', 'raw', 'dn', 'dn_biphasic')





