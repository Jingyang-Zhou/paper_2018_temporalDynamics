%dn_mkFigure4_2

%% LOAD DATA AND PARAMETER FILES

dataLoc    = fullfile(dn_ECoG_RootPath, 'data');
dt_fName   = 'dn_preprocessedData.mat';
prm_fName  = 'dn_params.mat';

a = load(fullfile(dataLoc, dt_fName));
b = load(fullfile(dataLoc, prm_fName));

%% EXTRACT PARAMETERS

w   = b.prm.ecog.dn_bi.weight;
ecc = a.dt.ecog.ecc;

prd = b.prm.ecog.dn_bi.prd;
dt  = a.dt.ecog.bbts_roi;

nrois = 3;

%% BIN WEIGHTS OVER ECCENTRICITIES

idx = {};

for iroi =  1: nrois
    idx.l{iroi} = (ecc{iroi} < 5)  &  (ecc{iroi} > .5);
    idx.m{iroi} = (ecc{iroi} >= 5) & (ecc{iroi} < 10);
    idx.h{iroi} = ecc{iroi} >= 10; 
end

%% PLOT PARAMETERS

figure (1), clf

col = [0.8*ones(1, 3); 0.5 * ones(1, 3); 0.1 * ones(1, 3)];

for k = 1 : nrois
   plot(1, mean(w{k}(idx.l{k})), 'o', 'markerfacecolor', col(k, :), 'markersize', 10, 'markeredgecolor', col(k, :)), hold on 
   plot(2, mean(w{k}(idx.m{k})), 'o', 'markerfacecolor', col(k, :), 'markersize', 10, 'markeredgecolor', col(k, :))
   plot(3, mean(w{k}(idx.h{k})), 'o', 'markerfacecolor', col(k, :), 'markersize', 10, 'markeredgecolor', col(k, :))
end
ylim([0.4, 1]), xlim([0.5, 3.5]), 
set(gca, 'xtick', 1 : 3, 'xticklabel', {'<5', '5-10', '>10'}, 'fontsize', 16)
xlabel('eccentricity'), ylabel('weight'), title('weight(w)')

%% PLOT EXAMPLES (Figure S5)

t = 0.001 : 0.001 : 1.204;

norm_max = @(x) x./max(x);

figure (2), clf
for k = 1 : nrois
   for k1 = 1 : length(w{k})
       subplot(3, 15, (k-1) * 15 + k1)
       plot(t, norm_max(dt{k}(:, k1)), 'k-', 'linewidth', 1.5), hold on
       plot(t, norm_max(prd{k}(k1, :)), 'r-', 'linewidth', 1.5), box off, 
       set(gca, 'xtick', '', 'ytick', '')
       title(sprintf('ecc = %0.2f', ecc{k}(k1)))
   end
end

subplot(3, 15, 1), ylabel('V1')
subplot(3, 15, 16), ylabel('V2')
subplot(3, 15, 31), ylabel('V3')
