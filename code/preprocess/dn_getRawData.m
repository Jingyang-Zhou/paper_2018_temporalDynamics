function [tsds, stimNames, imNames, goodLabels, goodChannels, idx] = dn_getRawData()

%% INITIALIZE DATA EXTRACTION

run('~/matlab/git/vistaproj_ECoG_pRF/ecogPRFpath.m')

% Data from which subject
subjID   = 'jt'; 
subjNum  = 17;
s        = ecogCreateSession(subjNum);

[dorsalIndex, dorsalLabels] = labelJtDorsalRoi;

% load all channels
channels     = ecogGet(s, 'good channels');
labels       = ecogGet(s, 'labels');
goodChannels = [channels, dorsalIndex];
goodLabels   = [labels, dorsalLabels];

%% GET RAW DATA PREP

s            = ecogSet(s, 'current stim', 'textures');
s            = ecogSet(s, 'current run', 1);

stimuli      = ecogGet(s, 'stimulusnumber', 'textures');
onsets       = ecogGet(s, 'onsets', 'textures', 1, 'frames');
onsets       = onsets(1:2:end);
epoch_length = median(diff(onsets));
image_file   = ecogGet(s, 'vista images', [], []);
tmp          = load(image_file);
% im           = tmp.stimulus.images{1};

imNames         = {'wn' 'pn' 'bn' 'gr8' 'gr16' 'gr32' 'gr64'};
nImgPerCategory = length(stimuli)/length(imNames);
stimNames       = ceil(stimuli / nImgPerCategory);

%% GET RAW DATA

for ii = 1:length(goodChannels)
    ts(:,ii)  = ecogGet(s, 'ts', 'textures', 1, goodChannels(ii), 'CAR');
    fprintf('.'); drawnow()
end
fprintf('\n');

%% DOWN-SAMPLE DATA

srate      = ecogGet(s, 'fs');
ts = double(ts);
for ii = 1 : length(goodChannels)
    tsd(:,ii)  = resample(ts(:,ii), 1000, 1031);
    tsds(:,ii) = resample(tsd(:,ii), 1000, 2960);
end

scale        = (1000/1031)*(1000/2960);
newsrate     = srate*scale;
onsets       = round(onsets.*scale);
epoch_length = median(diff(onsets));

%% EPOCH THE DATA/INDEX

frontpad = 200; % 200 ms

idx = zeros(length(onsets), epoch_length + frontpad);

for ii = 1:length(onsets)
    idx(ii, :) = onsets(ii) - frontpad + (0:epoch_length-1 + frontpad);
end

end
