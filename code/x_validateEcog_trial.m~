% x_validate_trialTypes

%% LOAD DATA

dataLoc    = fullfile('/Volumes/server/Projects/Temporal_integration/DN_2018_code_data/', 'data');
dt_fName   = 'dn_preprocessedData.mat';
prm_fName  = 'dn_params.mat';

a = load(fullfile(dataLoc, dt_fName));
b = load(fullfile(dataLoc, prm_fName));

%% PRE-DEFINED PARAMETERS

prm     = b.prm;
dn      = prm.ecog.dn;

irfType = 'uniphasic';
t       = prm.ecog.t;
stim    = prm.ecog.stim;

ecog = a.dt.ecog;

nrois = 4;

normmax = @(x)(x-mean(x(1 : 200)))./max(x - mean(x(1 : 200)));

%% extract broadband time course for each stimulus type

idx = []; dt = [];

for k = 1 : 3
    idx(k, :) = find(a.dt.ecog.stimNm == k);
end

for iroi = 1 : nrois
    [~,roiidx] = intersect(ecog.chans, ecog.elec_roi{iroi});
    tmp1 = ecog.bb{1}(idx(1, :), :, roiidx);
    tmp2 = ecog.bb{1}(idx(2, :), :, roiidx);
    tmp3 = ecog.bb{1}(idx(3, :), :, roiidx);
    tmp  = (tmp1 + tmp2 + tmp3)./3;
    dt{iroi} = mean(tmp, 3);
end

%%  VISUALIZE

figure (1), clf
for iroi = 1 : nrois
    subplot(1, 4, iroi), set(gca, 'colororder', copper(30)), hold on
    plot(dt{iroi}')
end

%% create training and test data

train = {}; ntrials = 30;

for iroi = 1 : nrois
    for k = 1 : ntrials
        idx = find([1 : ntrials]~=k);
        train{iroi}(:, k) = normmax(mean(dt{iroi}(idx, :)));
    end
end

%% GRID FIT TO EACH TRAINING DATA

for iroi = 1 : nrois
    [seed{iroi}, seedR2{iroi}] = dn_gridFit(train{iroi}', dn, stim, t, irfType);
end

%% FINE FIT 

for iroi = 1 : nrois
    thisseed = [seed{iroi}, zeros(size(seed{iroi}, 1), 1)];
    [param{iroi}, prd{iroi}, r2{iroi}, exitflg{iroi}] = dn_fineFit(train{iroi}', stim, t, dn, thisseed, irfType);
end

%% LEAVE-ONE-OUT CROSS-VALIDATION ACROSS ELECTRODES

xr2 = [];

computer2 = @(prd, dt) 1 - sum((prd - dt).^2)/sum((dt - mean(dt)).^2);

for iroi = 1 : nrois
    % compute correlation coefficient between prd and the original data
    for k = 1 : ntrials
        %norm_prd = prd{iroi}(k, :)./max(prd{iroi}(k, :));
        xr2{iroi}(k) = computer2(normmax(prd{iroi}(k, :)), normmax(dt{iroi}(k, :)));
    end
end

%% VISUALIZE XR2

figure (2), clf
for iroi = 1 : nrois
   plot(iroi, median(xr2{iroi}), 'ro', 'markerfacecolor', 'r', 'markersize', 15), hold on 
   plot(iroi, xr2{iroi}, 'ko', 'markerfacecolor', 'k', 'markeredgecolor', 'w'), 
end
box off, ylim([0, 1]), set(gca, 'fontsize', 16), xlim([0.5, 4.5]), set(gca, 'xtick', 1 : nrois), 
set(gca, 'xticklabel', {'V1', 'V2', 'V3', 'anterior'}), title('x-validated R2')

%% VISUALIZE FITS

figure (3), clf
for iroi = 1 : nrois
    for k = 1 : ntrials
        subplot(4, ntrials, (iroi-1)*ntrials + k), 
        plot(t, normmax(dt{iroi}(k, :)), 'k-', 'linewidth', 2), hold on
        plot(t, normmax(prd{iroi}(k, :)), 'r-', 'linewidth', 2), axis tight, box off
        set(gca, 'xtick', [0.2, 0.7], 'ytick', [0, 1]), axis off
    end
end

%% save

prm = b.prm;
prm.ecog.dn.xr2.trial_r2  = xr2;
prm.ecog.dn.xr2.trial_prd = prd;
prm.ecog.dn.xr2.trial_dt  = dt;

save(fullfile(dataLoc, prm_fName), 'prm')