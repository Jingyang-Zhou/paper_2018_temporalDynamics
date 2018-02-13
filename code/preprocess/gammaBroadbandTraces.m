function [bbmatrix, stimNames, imNames, goodLabels, goodChannels, srate] = gammaBroadbandTraces(subjID, channels, bbmethod)

%%

% INPUTS: 

% channels : either empty (load all channels), or a vector, each
%           entry represents an electrode.

% subjID   : either "jt" or "rw" for now

% OUTPUTS: 
% subjID = 'jt';


%%
% load subject channels and labels
figureOn = 0;

run('~/matlab/git/vistaproj_ECoG_pRF/ecogPRFpath.m')

if ~exist('subjID', 'var'), 
    subjID = 'jt'; subjNum = 17;
elseif strcmpi(subjID, 'jt'),
    subjNum = 17; 
end

s                           = ecogCreateSession(subjNum);
[dorsalIndex, dorsalLabels] = labelJtDorsalRoi;

% lead load all channels or selected channels

if ~exist('channels', 'var') || isempty(channels)
    channels              = ecogGet(s, 'good channels');
    labels                = ecogGet(s, 'labels');
    goodChannels          = [channels, dorsalIndex];
    goodLabels            = [labels, dorsalLabels];
    
else
    allChannels           = s.electrodes;
    visualChannels        = [ecogGet(s, 'good channels'), dorsalIndex];
    labels                = [ecogGet(s, 'labels'), dorsalLabels];
    
    assert(ismember(channels, allChannels), 'Channel input is not on the electrode list');
    
    if all(ismember(channels, visualChannels))
        [~, index]            = ismember(channels, visualChannels);
        goodLabels            = allLabels(1, index);
    else
        goodLabels            = channels;
    end
    goodChannels          = channels;
end

synthetic_data = false;


%% get raw data prep

s            = ecogSet(s, 'current stim', 'textures');
s            = ecogSet(s, 'current run', 1);

stimuli      = ecogGet(s, 'stimulusnumber', 'textures');
onsets       = ecogGet(s, 'onsets', 'textures', 1, 'frames');
onsets       = onsets(1:2:end);
epoch_length = median(diff(onsets));
image_file   = ecogGet(s, 'vista images', [], []);
tmp          = load(image_file);
% im           = tmp.stimulus.images{1};

imNames   = {'wn' 'pn' 'bn' 'gr8' 'gr16' 'gr32' 'gr64'};
nImgPerCategory = length(stimuli)/length(imNames);
stimNames   = ceil(stimuli / nImgPerCategory);

%% Get the raw data 

for ii = 1:length(goodChannels)
    texture_ts(:,ii)  = ecogGet(s, 'ts', 'textures', 1, goodChannels(ii), 'CAR');
    fprintf('.'); drawnow()
end
fprintf('\n');

%% down-sample data 

srate      = ecogGet(s, 'fs');
texture_ts = double(texture_ts);
for ii = 1 : length(goodChannels)
    texture_tsd(:,ii)  = resample(texture_ts(:,ii), 1000, 1031);
    texture_tsds(:,ii) = resample(texture_tsd(:,ii), 1000, 2960);
end
scale    = (1000/1031)*(1000/2960);
newsrate = srate*scale;
onsets   = round(onsets.*scale);
epoch_length = median(diff(onsets));

%% For simulated rather than real data 

if synthetic_data
    simulated_baseline = randn(size(texture_tsds));
    simulated_signal   = randn(size(texture_tsds));
    irf_concatenated   = zeros(size(simulated_baseline));
    
    tau1 = 50; % ms
    onset_delay = 200; % ms
    
    % IRF
    t = 0:999;
    irf = (t-onset_delay) .* exp(-(t-onset_delay)/tau1);
    irf(t<onset_delay) = 0;
    irf = irf/sum(irf) * newsrate;
    irf = repmat(irf, [size(simulated_baseline,2), 1])';
    % Concatenated IRF
    for ii = 1:length(onsets)
        irf_concatenated((1:length(t))+onsets(ii),:) = irf;
    end
    
    simulated_signal = simulated_baseline + ...
        simulated_signal .* irf_concatenated;
    
    texture_tsds = simulated_signal;
end

%% Extract broadband signal 

%   define bands for subsequent band pass filtering
band_rg  = [80 200];
band_w   = 10;
lb       = band_rg(1):band_w:band_rg(2)-band_w;
ub       = lb+band_w;
bands   = [lb; ub]';
disp(bands)


%    filter
whiten = @(x) bsxfun(@plus, zscore(x), mean(x));
for c = 1:size(texture_tsds,2)
    fprintf('.'); drawnow()
    signal.raw = double(texture_tsds(:,c));
    signal.bp_multi  = zeros(length(signal.raw),size(bands,1));
    for ii = 1:size(bands,1)
        [signal.bp_multi(:,ii)]=butterpass_eeglabdata(signal.raw,bands(ii,:),newsrate);
    end
    
    % THIS PART MAY CHANGE IN THE FUTURE **************************
    
    signal.bb = sum(whiten((abs(hilbert(signal.bp_multi)))),2);
    
    texture_bb(:,c) = signal.bb;
end
fprintf('\n');

% Epoch the data  -------------------------------------------------------

frontpad = 200; % 200 ms

texture_tsmatrix = zeros(epoch_length + frontpad, length(onsets), length(goodChannels));
bbmatrix = zeros(epoch_length + frontpad, length(onsets), length(goodChannels));

for ii = 1:length(onsets)
    %idx = onsets(ii)+(0:epoch_length-1);
    
    idx = onsets(ii) - frontpad + (0:epoch_length-1 + frontpad);
    texture_tsmatrix(:,ii,:) = texture_tsds(idx,:);
    bbmatrix(:,ii,:) = texture_bb(idx,:);
end

%% plot and save channels

if figureOn
    fH = figure; pos = get(fH, 'Position'); pos([3 4]) = [750 1200];
    set(fH, 'Color', 'w', 'Position', pos);
    
    for ii = 1:length(goodChannels)
        t = (1:epoch_length + frontpad)/newsrate;
        y1 = mean(texture_tsmatrix(:,:,ii),2);
        y1 = y1 / max(abs(y1));
        y2 = mean(bbmatrix(:,:,ii),2);
        y2 = y2 - min(y2);
        y2 = y2 / max(y2);
        plot(t, y1, t,y2)
        title(sprintf('Electrode %d\t%s %d', goodChannels(ii), goodLabels{ii}, ii))
        set(gca, 'XTick', 0:.1:1, 'YLim', [-1 1], 'XGrid', 'on')
        %fname = fullfile(temporalRootPath, 'scratch', sprintf('subj%02d_ECoG_Channel_%03d.eps', subjNum, goodChannels(ii)));
        %hgexport(fH, fname);
        waitforbuttonpress
    end
end

%% select channels

% finaltexture_bbmatrix = [];
% finalGoodLabels       = {};
% 
% on  = 201 : 700;
% off = [1 :200, 700 : 1200];
% k   = 0;
% idx = [];
% 
% for ii = 1 : length(goodChannels)
%     tmp = mean(texture_bbmatrix(:, :, ii), 2);
%     if std(tmp(on)) > 2*std(tmp(off)), 
%         k = k + 1;
%         %finaltexture_bbmatrix(:, :, k) = texture_bbmatrix(:, :, ii);
%         finalGoodLabels{k}             = goodLabels{ii};
%         idx(k) = ii;
%     else
%         disp(ii)
%     end
% end

%%
end

