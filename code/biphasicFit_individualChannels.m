% biphasic fit to individual electrodes

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

irfType = 'biphasic';
t       = prm.ecog.t;
stim    = prm.ecog.stim;
nrois   = length(dt) - 1;

ecc  = a.dt.ecog.ecc;

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

%% Biphasic fit to data

seedval  = [0.05, 0.6, 0.1, 0.05, 0];
bound.lb = [0.01, 0, 0.0001,  0.01, 0];
bound.ub = [1, 1, 1,  .5, 0.5];

nrois = 3;


for iroi = 1 : nrois
    seed = repmat(seedval, [size(dt{iroi}, 2), 1]);
    [param{iroi}, prd{iroi}, r2{iroi}] = ...
        dn_fineFit(dt{iroi}', stim, t, bound, seed, irfType);
end

% VISUALIZE THE WEIGHT

w = [];
for k = 1 : nrois
    w{k} = param{k}(:, 2);
end

% BIN WEIGHT OVER ECCENTRICITIES

idx = {};

for iroi =  1: nrois
    idx.l{iroi} = (ecc{iroi} < 5)  &  (ecc{iroi} > .5);
    idx.m{iroi} = (ecc{iroi} >= 5) & (ecc{iroi} < 10);
    idx.h{iroi} = ecc{iroi} >= 10; 
end

% PLOT PARAMETERS

figure (1), clf

for k = 1 : nrois
   plot(1, mean(w{k}(idx.l{k})), 'ko', 'markerfacecolor', 'k'), hold on 
   plot(2, mean(w{k}(idx.m{k})), 'bo', 'markerfacecolor', 'b')
   plot(3, mean(w{k}(idx.h{k})), 'ro', 'markerfacecolor', 'r')
end
%ylim([0.4, 1])

norm_max = @(x) x./max(x);

figure (2), clf
for k = 1 : nrois
    for k1 = 1 : length(w{k})
        subplot(3, 15, (k-1) * 15 + k1)
       plot(t, norm_max(dt{k}(:, k1))), hold on
       plot(t, norm_max(prd{k}(k1, :)))
    end
end
