% dn_fitToECoGFrequencyBands

%% load data

dataLoc = fullfile(temporalRootPath, 'data');
fName   = 'dn_data.mat';

a = load(fullfile(dataLoc, fName));

stimTypes = {'wn' 'pn' 'bn' 'gr8' 'gr16' 'gr32' 'gr64'};

stimToUse = [1 : 3, 4, 5, 6, 7];
ts_raw    = a.raw.tsmatrix;
type_idx  = a.raw.stimuliNames;

% locate V1 labels
v3Labels = []; idx = 0;

for k = 1 : length(a.raw.goodLabels)
    if contains(a.raw.goodLabels{k}, 'V1')
        idx = idx + 1;
        v3Labels(idx) = k;
    end
end

%% useful function 

normMax_range = @(x, range) (x - mean(x(range)))./max(x - mean(x(range)));

%% extract broadband time course w.r.t. each stimulus types: prep

ts_types = {};
bb_types = {};

% extract different frequency bands
f_Rng  = [1, 200]; bands_indi = {}; idx = 0; srate = 1000;

for k = f_Rng(1) : 30 : f_Rng(2)
    idx = idx + 1; 
    bands_indi{idx} = {[k, k + 30], 10};
    bands_indi_str{idx} = sprintf('[%d, %d]', k, k + 30);
end

%% extract broadband

ts_types = {};

for k = stimToUse
    ts_types{k } = ts_raw(:, find(type_idx == k), :);
    % extract broadband
    for k1 = 1 : length(bands_indi)
        for k2 = 1 : size(ts_types{1}, 3) % number of electrodes
            bb_types{k, k1}(:, :, k2) = normMax_range(extractBroadband(squeeze(ts_types{k}(:, :, k2)), srate, 4, bands_indi{k1}), [1 : 200]); 
        end
    end
    disp(k)
end
%% average V1 broadband response

bb_types_v1 = [];

for k = 1 : length(stimToUse)
    for k1 = 1 : length(bands_indi)
        tmp = median(median(bb_types{k, k1}(100 : 1000, :, v3Labels), 2), 3);
        bb_types_v1(k, k1, :) = normMax_range(tmp, [1 : 100]);
    end
end

%% visualize broadband response to each grating stimulus

figure (1), clf
for k = 1 : length(stimToUse)
   subplot(7, 1, k), cla, set(gca, 'ColorOrder', copper(7)), hold on
   plot(squeeze(bb_types_v1(k, 2 : end, :))'), axis tight, 
   title(stimTypes{k}), ylim([-0.2, 1])
   %legend(bands_indi_str{2 : end}), set(gca, 'fontsize', 12)
end