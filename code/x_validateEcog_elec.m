% do cross-one electrode out cross validation

% I'm planning to cross-validate across electrodes, and probably I'm going
% to fit an additional scaling factor to the predicted electrodes.

% First fit uniphasic, then fit biphasic;

%% LOAD DATA

dataLoc    = fullfile('/Volumes/server/Projects/Temporal_integration/DN_2018_code_data/', 'data');
dt_fName   = 'dn_preprocessedData.mat';
prm_fName  = 'dn_params.mat';

a = load(fullfile(dataLoc, dt_fName));
b = load(fullfile(dataLoc, prm_fName));

%% PRE-DEFINED PARAMETERS

dt      = a.dt.ecog.bbts_roi;
prm     = b.prm;
dn      = prm.ecog.dn;

irfType = 'uniphasic';
t       = prm.ecog.t;
stim    = prm.ecog.stim;
nrois   = length(dt);

%% NORMALIZE EACH TIME COURSE UP TO 1

dt_norm = {};

for iroi = 1 : nrois
    for k = 1 : size(dt{iroi}, 2)
        dt_norm{iroi}(:, k) = dt{iroi}(:, k)./max(dt{iroi}(:, k));
    end
end

%% SPECULATE DATA

figure (1), clf
for iroi = 1 : nrois
    subplot(2, 4, iroi), set(gca, 'colororder', copper(size(dt{iroi}, 2))), hold on
    plot(dt{iroi}, 'linewidth', 2), axis tight
    
    subplot(2, 4, iroi + 4), set(gca, 'colororder', copper(size(dt{iroi}, 2))), hold on
    plot(dt_norm{iroi}, 'linewidth', 2), axis tight
end

%% CREATE TRAINING AND TEST DATA

train = {};

for iroi = 1 : nrois
    nelec = size(dt_norm{iroi}, 2);
    for k = 1 : nelec
        idx = find([1 : nelec]~=k);
        train{iroi}(:, k) = mean(dt_norm{iroi}(:, idx), 2);
    end
end
% train model on the training data, and test the model on the left out
% (original) data.

%% GRID FIT AND FINE FIT TO EACH TRAINING DATA

for iroi = 1 : nrois
    [seed{iroi}, seedR2{iroi}] = dn_gridFit(train{iroi}', dn, stim, t, irfType);
end

%% FINE FIT TO THE TRAINING DATA

for iroi = 1 : nrois
    thisseed = [seed{iroi}, zeros(size(seed{iroi}, 1), 1)];
    [param{iroi}, prd{iroi}, r2{iroi}, exitflg{iroi}] = dn_fineFit(train{iroi}', stim, t, dn, thisseed, irfType);
end

%% LEAVE-ONE-OUT CROSS-VALIDATION ACROSS ELECTRODES

xr2 = [];

computer2 = @(prd, dt) 1 - sum((prd - dt).^2)/sum((dt - mean(dt)).^2);

for iroi = 1 : nrois
    nelec = size(dt_norm{iroi}, 2);
    % compute correlation coefficient between prd and the original data
    for k = 1 : nelec
        norm_prd = prd{iroi}(k, :)./max(prd{iroi}(k, :));
        xr2{iroi}(k) = computer2(norm_prd, dt_norm{iroi}(:, k)');
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

normmax = @(x)x./max(x);

figure (3), clf
for iroi = 1 : nrois
    for k = 1 : size(dt_norm{iroi}, 2)
        subplot(4, 15, (iroi-1)*15 + k), 
        plot(t, dt_norm{iroi}(:, k), 'k-', 'linewidth', 2), hold on
        plot(t, normmax(prd{iroi}(k, :)), 'r-', 'linewidth', 2), axis tight, box off
        set(gca, 'xtick', [0.2, 0.7], 'ytick', [0, 1]), axis off
    end
end

%% save

prm = b.prm;
prm.ecog.dn.xr2.elec_r2  = xr2;
prm.ecog.dn.xr2.elec_prd = prd;

save(fullfile(dataLoc, prm_fName), 'prm')