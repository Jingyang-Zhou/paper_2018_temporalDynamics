% dn_fitDNECoG2fMRI

%% SAVE FIGURE KNOB

saveFigure = 0;

%% LOAD fMRI AND ECOG DATA

dataLoc = fullfile('/Volumes/server/Projects/Temporal_integration/DN_2018_code_data/', 'data');
fName   = 'trf_fmriData.mat';
%fName1  = 'dn_params.mat';
fName2   = 'dn_data.mat';

% LOAD fMRI DATA
a = load(fullfile(dataLoc, fName));
fmriData = squeeze(mean(a.fmri1.data));

% LOAD ECOG DATA
b = load(fullfile(dataLoc, fName2));

% Make fMRI stimulus
[stimulus, time] = mk_fMRIStimulus('with0');

%% PRE-DEFINED / DERIVED PARAMETERS

% fMRI PARAMETERS ---------------------------------------------------------
nBoots    = 100;
nfmriRois = 9;

fmri_2fit(1 : 3, :, :) = fmriData(1 : 3, :, :); % V1-V3 data

tmp = fmriData(4, :, :);

for k = 5 : nfmriRois
    tmp  = tmp + fmriData(k, :, :);
end
fmri_2fit(4, :, :) = squeeze(tmp)./6; % 6 ROIs get combined into 1

% ECOG PARAMETERS ---------------------------------------------------------

%ecogPrm = b.prm.ecog.dn.param;

ecogPrm = b.dn_biphasic.biprm;

roiNm = {'V1', 'V2', 'V3'};
nrois = length(roiNm);

for iroi = 1 : length(roiNm)
    m_ecogPrm(iroi, :) = median(ecogPrm{iroi});
end

%% COMPUTE MEAN AND STD OF THE FMRI DATA

mfmri_2fit = []; sfmri_2fit = [];

for iroi = 1 : nrois
    mfmri_2fit(iroi, :) = squeeze(median(fmri_2fit(iroi, :, :), 2));
    sfmri_2fit(iroi, :, :) = squeeze(prctile(fmri_2fit(iroi, :, :), [25, 75], 2));
end

%% GENERATE DN PREDICTIONS TO ALL FMRI STIMULI

scale_max = @(x) x./max(x(:));

dn_pred = [];

for iroi = 1 : nrois
    prm = [m_ecogPrm(iroi, 1 : 3), 2, m_ecogPrm(iroi, 4 : 5), 1];
    %prm(3) = 0.0001;
    dn_pred(iroi, :, :) = scale_max(dn_DNmodel(prm, stimulus, time));
end

%% FIT ECOG TO FMRI

scale = [];

fitType1 = 'linear';
%fitType2 = 'sqrt';

seed = 0.001;

linrsp = [1, 2, 4, 8, 16, 32, 16, 16, 16, 16, 16, 16, 0];

for iroi = 1 : nrois
    ecogrsp(iroi, :) = sum(dn_pred(iroi, :, :), 3);
    mrirsp  = mfmri_2fit(iroi, :);
    scale(1, iroi) = fminsearch(@(x) dn_fitEcog2fMRI(x, ecogrsp(iroi, :), mrirsp, fitType1), seed);
    %scale(2, iroi) = fminsearch(@(x) dn_fitEcog2fMRI(x, ecogrsp(iroi, :), mrirsp, fitType2), seed);
    
    % fit the linear model and its scale to the fMRI data
    linscale(iroi) = fminsearch(@(x) dn_fitEcog2fMRI(x, linrsp, mrirsp, fitType1), seed);
    % compute ECoG prediction
    ecogPrd1(iroi, :) = scale(1, iroi) * ecogrsp(iroi, :);
    
    % compute r2
    r2(iroi, 1) = corr(ecogPrd1(iroi, :)', mrirsp').^2;
    
    % compute the linear r2
    linr2(iroi) = corr(linscale(iroi) * linrsp', mrirsp').^2;
end


%% VISUALIZE fMRI DATA

x1 = [0, 1, 2, 4, 8, 16, 32];

fg1 = figure (1); clf

iroi = 1;

for jj = 1:2
    switch jj
        case 1, mriOrder = [13, 1 :  6];
        case 2, mriOrder = [5 7:12];
    end
    
    subplot(1, 2,  jj)
    plot(x1, mfmri_2fit(iroi, mriOrder), 'ko', 'markersize', 9,  'markerfacecolor', 'k'), hold on
    for k = 1 : length(mriOrder)
        low = sfmri_2fit(iroi, 1, mriOrder(k));
        high = sfmri_2fit(iroi, 2, mriOrder(k));
        plot(x1(k)*[1 1], [low, high], '-', 'linewidth', 3, 'color', 0.8 * ones(1, 3))
    end
    
    ylim([-.05, 0.6]), xlim([-2, 34]), box off, set(gca, 'xtick', [0, 4, 16, 32], 'xticklabel', [0, 67, 267, 533])
    
    title(roiNm{iroi}), set(gca, 'fontsize', 14)
    
    plot(x1, linscale(iroi).*linrsp(mriOrder), 'b-', 'linewidth', 3)
    plot(x1, ecogPrd1(iroi, mriOrder), 'r-', 'linewidth', 3)
    set(gca, 'xaxislocation', 'origin'), xlabel('time (ms)'), ylabel('% BOLD')
end


