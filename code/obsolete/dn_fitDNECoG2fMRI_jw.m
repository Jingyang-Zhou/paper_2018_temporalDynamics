% dn_fitDNECoG2fMRI

%% SAVE FIGURE KNOB

saveFigure = 0;

%% LOAD fMRI AND ECOG DATA

dataLoc = fullfile('/Volumes/server/Projects/Temporal_integration/DN_2018_code_data/', 'data');
fName   = 'trf_fmriData.mat';
fName1  = 'dn_params.mat';

fName2 = 'dn_data.mat';

% LOAD fMRI DATA
a = load(fullfile(dataLoc, fName));
fmriData = squeeze(mean(a.fmri1.data));

% LOAD ECOG DATA
b = load(fullfile(dataLoc, fName1));
c = load(fullfile(dataLoc, fName2));

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
nrois   = 4;
ecogPrm = b.prm.ecog.dn.param;

roiNm = b.prm.ecog.roiNm;

% combine biphasic ecog parameters
biprm = [];
biprm = c.dn_biphasic.biprm;
for k = 1 : 3
    m_biprm(k, :) = median(biprm{k});
end

%% COMPUTE MEAN AND STD OF THE FMRI DATA

mfmri_2fit = []; sfmri_2fit = [];

for iroi = 1 : nrois
    mfmri_2fit(iroi, :) = squeeze(median(fmri_2fit(iroi, :, :), 2));
    sfmri_2fit(iroi, :, :) = squeeze(prctile(fmri_2fit(iroi, :, :), [25, 75], 2));
end

%% COMPUTE THE MEDIAN OF THE ECOG PARAMETERS FOR EACH ROI

mprm = [];

for iroi = 1 : nrois
    idx = (iroi - 1) * nBoots + 1 : iroi * nBoots;
    mprm(iroi, :) = median(ecogPrm(idx, :));
end

%% GENERATE DN PREDICTIONS TO ALL FMRI STIMULI

scale_max = @(x) x./max(x(:));

dn_pred = [];

nrois = 3;

for iroi = 1 : nrois
    prm = [mprm(iroi, 1), 0, mprm(iroi, 2), mprm(iroi, 3), mprm(iroi, 4 : 5), 1];
    prm = [m_biprm(iroi,1 : 3), 2, m_biprm(iroi, 4 : 5), 1];
    % prm(3) = 0.0001;
    dn_pred(iroi, :, :) = scale_max(dn_DNmodel(prm, stimulus, time));
end

%% FIT ECOG TO FMRI

scale = [];

fitType1 = 'linear';
fitType2 = 'sqrt';

seed = 0.001;

linrsp = [1, 2, 4, 8, 16, 32, 16, 16, 16, 16, 16, 16, 0];

for iroi = 1 : nrois
    ecogrsp(iroi, :) = sum(dn_pred(iroi, :, :), 3);
    mrirsp  = mfmri_2fit(iroi, :);
    scale(1, iroi) = fminsearch(@(x) dn_fitEcog2fMRI(x, ecogrsp(iroi, :), mrirsp, fitType1), seed);
    scale(2, iroi) = fminsearch(@(x) dn_fitEcog2fMRI(x, ecogrsp(iroi, :), mrirsp, fitType2), seed);
    
    % fit the linear model and its scale to the fMRI data
    linscale(iroi) = fminsearch(@(x) dn_fitEcog2fMRI(x, linrsp, mrirsp, fitType1), seed);
    % compute ECoG prediction
    ecogPrd1(iroi, :) = scale(1, iroi) * ecogrsp(iroi, :);
    ecogPrd2(iroi, :) = scale(2, iroi) * sqrt(ecogrsp(iroi, :));
    % compute r2
    r2(iroi, 1) = corr(ecogPrd1(iroi, :)', mrirsp').^2;
    r2(iroi, 2) = corr(ecogPrd2(iroi, :)', mrirsp').^2;
    % compute the linear r2
    linr2(iroi) = corr(linscale(iroi) * linrsp', mrirsp').^2;
end

% linscale, linrsp
%% COMPUTE TEMPORAL SUMMAITON RATIO

% mfmri_2fit

% ecogPrd1

% ecogPrd2

sumRatio = [];

for iroi = 1 : nrois
    for k = 1 : 5
        sumRatio(1, iroi, k) = mfmri_2fit(iroi, k + 1) / (mfmri_2fit(iroi, k) * 2);
        sumRatio(2, iroi, k) = ecogPrd1(iroi, k + 1) / (ecogPrd1(iroi, k) * 2);
        sumRatio(3, iroi, k) = ecogPrd2(iroi, k + 1) / (ecogPrd2(iroi, k) * 2);
    end
end

msumRatio = squeeze(median(sumRatio, 3));


%% PLOT SUMMATION RATIO

figure (2), clf, set(gca, 'colororder', copper(3)), hold on
%plot(msumRatio', 'o')
plot(msumRatio(1, :), 'ro', 'markersize', 7, 'markerfacecolor', 'r'), hold on
plot(msumRatio(2, :)', 'bo', 'markersize', 7, 'markerfacecolor', 'b')
plot(msumRatio(3, :)', 'ko', 'markersize', 7, 'markerfacecolor', 'k')
ylim([0.5, 1.2]), xlim([0.5, nrois + 0.5]),
set(gca, 'xtick', 1 : nrois, 'xticklabel', roiNm, 'fontsize', 14)
legend('data summation ratio', 'linear prediction', 'sqrt DN predicted summation ratio'),


%% VISUALIZE PREDICTIONS

% figure (100), clf, set(gca, 'colororder', copper(2)), hold on
% for k = 1 : 13
%     plot(time, squeeze(dn_pred(3, [5, 12], :)), 'linewidth', 3),
% end
% xlim([0, 1.5]), legend('rsp to 267ms single pulse', 'rsp. to two 134ms pulses with 533ms gap'),
% set(gca, 'fontsize', 14)

%% VISUALIZE fMRI DATA

%mriOrder = [13, 1 :  6, 5, 7 : 12];


x1 = [0, 1, 2, 4, 8, 16, 32];

fg1 = figure (1); clf
for iroi = 1 : nrois
    
    for jj = 1:2
        switch jj 
            case 1, mriOrder = [13, 1 :  6];
            case 2, mriOrder = [5 7:12];
        end
        
    subplot(2, 4, iroi + 4*(jj-1))
    
    plot(x1, mfmri_2fit(iroi, mriOrder), 'ko', 'markersize', 7, 'linewidth', 1), hold on
    
    for k = 1 : length(mriOrder)
        low = sfmri_2fit(iroi, 1, mriOrder(k));
        high = sfmri_2fit(iroi, 2, mriOrder(k));
        plot(x1(k)*[1 1], [low, high], 'k-')
    end
    %xlim([0.5, length(mriOrder) + 0.5]),
    ylim([-.05, 0.6]), box off
    %xlim(7+[0.5, 7 + 0.5]), ylim([-.05, 0.6]), box off
    
%     set(gca, 'xaxislocation', 'origin'), set(gca, 'xtick', [1 : length(mriOrder)], 'xticklabel', '')
%     text(7, 0.06, sprintf('r^2: %.2f', r2(iroi, 1)))
%     text(11, 0.06, sprintf('r^2: %.2f', r2(iroi, 2)), 'color', 'r')
    title(roiNm{iroi}), set(gca, 'fontsize', 14)
    
    
    plot(x1, linscale(iroi).*linrsp(mriOrder), 'g-', 'linewidth', 2)
    plot(x1, ecogPrd1(iroi, mriOrder), '-.', 'linewidth', 2, 'color', 0 * ones(1, 3))
    
    % PLOT LINEAR NEURONAL PREDICTION
    %plot(x1, linscale(iroi).*linrsp(mriOrder), 'g-', 'linewidth', 2),
    
    % PLOT LINEAR ECOG PREDICTIONS
    %plot(x1, ecogPrd1(iroi, mriOrder(1 : 7)), '-.', 'linewidth', 2, 'color', 0 * ones(1, 3)), hold on
    
    % PLOT SQRT ECOG PREDICTIONS
    
    % plot(x1, ecogPrd2(iroi, mriOrder), 'k-', 'linewidth', 2),
    % plot(x2, ecogPrd2(iroi, mriOrder(8 : 14)),'k-',  'linewidth', 2),
    end
end

fg1.Position = [0, 100, 1000, 600];

%% VISUALIZE fMRI DATA 2

mriOrder = [13, 1 :  6, 5, 7 : 12];

fg2 = figure (99); clf
for iroi = 1 : nrois
    subplot(1, 4, iroi)
    plot(mfmri_2fit(iroi, mriOrder), 'ko', 'markersize', 7, 'linewidth', 1), hold on
    for k = 1 : length(mriOrder)
        low = sfmri_2fit(iroi, 1, mriOrder(k));
        high = sfmri_2fit(iroi, 2, mriOrder(k));
        plot([k, k], [low, high], 'k-')
    end
    % PLOT LINEAR NEURONAL PREDICTION
    plot(1 : 7, linscale(iroi).*linrsp(mriOrder(1 : 7)), 'g-', 'linewidth', 2),
    plot(8 : 14, linscale(iroi).*linrsp(mriOrder(8 : 14)), 'g-', 'linewidth', 2)
    % PLOT LINEAR ECOG PREDICTIONS
%     plot(1 : 7, ecogPrd1(iroi, mriOrder(1 : 7)), '-.', 'linewidth', 2, 'color', 0 * ones(1, 3)), hold on
%     plot(8 : 14, ecogPrd1(iroi, mriOrder(8 : 14)), '-.', 'linewidth', 2, 'color', 0 * ones(1, 3))
%     % PLOT SQRT ECOG PREDICTIONS
%     plot(1 : 7, ecogPrd2(iroi, mriOrder(1 : 7)), 'k-', 'linewidth', 2),
%     plot(8 : 14, ecogPrd2(iroi, mriOrder(8 : 14)),'k-',  'linewidth', 2),
    xlim([0.5, length(mriOrder) + 0.5]), ylim([-.05, 0.45]), box off
    set(gca, 'xaxislocation', 'origin'), set(gca, 'xtick', [1 : length(mriOrder)], 'xticklabel', '')
%     text(7, 0.06, sprintf('r^2: %.2f', r2(iroi, 1)))
%     text(11, 0.06, sprintf('r^2: %.2f', r2(iroi, 2)), 'color', 'r')
    title(roiNm{iroi}), set(gca, 'fontsize', 14)
end

fg2.Position = [0, 100, 1000, 100];

%% SAVE FIGURE

if saveFigure,
    fgNm = 'fg_ecogDNFit2MRI_1';
    figLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
    printnice(1, 0, figLoc, fgNm);
    
    fgNm = 'fg_ecogDNFit2MRI_2';
    figLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
    printnice(2, 0, figLoc, fgNm);
end
