% dn_mkFigureS3

%% LOAD DATA AND PARAM FILES

dataLoc    = fullfile(dn_ECoG_RootPath, 'data');
dt_fName   = 'dn_preprocessedData.mat';
prm_fName  = 'dn_params.mat';

a = load(fullfile(dataLoc, dt_fName));
b = load(fullfile(dataLoc, prm_fName));

%% VISUALIZE FITS (ACROSS TRIALS)

normmax = @(x)(x - mean(x(1 : 200)))./max(x-mean(x(1 : 200)));

trial_dt  = b.prm.ecog.dn.xr2.trial_dt;
trial_prd = b.prm.ecog.dn.xr2.trial_prd;
nrois     = 4;

t = 0.001 : 0.001 : 1.204;

figure (1), clf
for iroi = 1 : nrois
    for k = 1 : size(trial_dt{iroi}, 1)
        subplot(4, 30, (iroi-1)*30 + k), 
        plot(t, normmax(trial_dt{iroi}(k, :)), 'k-', 'linewidth', 2), hold on
        plot(t, normmax(trial_prd{iroi}(k, :)), 'r-', 'linewidth', 2), axis tight, box off
        set(gca, 'xtick', [0.2, 0.7], 'ytick', [0, 1]), axis off
    end
end

%% VISUALIZE FITS (ACROSS ELECTRODES)

elec_dt  = a.dt.ecog.bbts_roi;
elec_prd = b.prm.ecog.dn.xr2.elec_prd;

figure (2), clf
for iroi = 1 : nrois
    for k = 1 : size(elec_prd{iroi}, 1)
        subplot(4, 15, (iroi-1)*15 + k), 
        plot(t, normmax(elec_dt{iroi}(:, k)), 'k-', 'linewidth', 2), hold on
        plot(t, normmax(elec_prd{iroi}(k, :)), 'r-', 'linewidth', 2), axis tight, box off
        set(gca, 'xtick', [0.2, 0.7], 'ytick', [0, 1]), axis off
    end
end