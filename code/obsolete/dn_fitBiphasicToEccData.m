% dn_fitBiphasicToEccData

% DEPENDENCIES: ecogFitPRFPix2Degl

% modification to params
% For biphasic fine fit, the lower bounds and upper bounds are changed
% for  biphasic coarse fit, tau1steps, tau2steps and weightsteps are
% changed

%% load data

dataLoc = fullfile(temporalRootPath, 'data');
fName   = 'dn_data.mat';

a       = load(fullfile(dataLoc, fName));
param   = a.param;

nrois = length(param.roiNm);
roiTS = param.roiTS;

%% check electrode eccentricities

b = load('/Volumes/server/Projects/ECoG/Figures/pRF_data/Power/oneBar/CAR/Subj17_exp_bb.mat');
x0 = b.params.params(:, 1);
y0 = b.params.params(:, 2);
s0 = b.params.params(:, 3);

subj = 17;

for iroi = 1 : nrois
    idx = param.elecIdxRelAll{iroi};
    [x{iroi}, y{iroi}, s{iroi}] = ecogFitPRFPix2Deg(subj, x0(idx), y0(idx), s0(idx));
    ecc{iroi} = sqrt(x{iroi}.^2 + y{iroi}.^2);
end

param.elecEcc = ecc;
%% creat bins

bin(1, :) = [0, 5];
bin(2, :) = [5, 10];
bin(3, :) = [10, inf];

%% assign eccntricity bins

eccBinIdx = [];

for iroi = 1 : nrois
    for k = 1 : size(roiTS{iroi}, 2)
        if  ecc{iroi}(k) < 5
            eccBinIdx{iroi}(k) = 1;
        elseif (ecc{iroi}(k) >= 5) & (ecc{iroi}(k) < 10)
            eccBinIdx{iroi}(k) = 2;
        elseif ecc{iroi}(k) >= 10
            eccBinIdx{iroi}(k) = 3;
        end
    end
end

%% compute offset response index - coarse fit

% coarse fit monophasic model to the first 700 ms of the time courses
t = linspace(0.001, 0.7, 700);
irfType = 'uniphasic';

stim = zeros(1, 700); 
stim(200 : 700) = 1;

for iroi = 1 : 3%nrois
    iroi
   for k = 1 : size(roiTS{iroi}, 2)
       k
       data = roiTS{iroi}(1 : 700, k)';
       [modelSeed{iroi}(k, :), seedR2{iroi}(k)] = ...
           dn_gridFit(data, param, stim, t, irfType);
   end
end

%% compute offset response index - fine fit

% if skip the grid fit ----------------------------
bi = a.dn_biphasic;
modelSeed = bi.ORImodelSeed;
irfType = 'uniphasic';
roiTS = a.param.roiTS';

t1 = 0.001 : 0.001 : 1.20;
stim1 = zeros(1, 1200); stim1(200 : 700) = 1;

for iroi = 1 : 3 %nrois
    iroi
    for k = 1 : size(roiTS{iroi}, 2)
        k
        % seed = [modelSeed{iroi}(k, 1),0.9, modelSeed{iroi}(3 : 4), 0];
        seed = [modelSeed{iroi}(k, :), 0];
        data = roiTS{iroi}(1 : 1200, k)';
        [prm{iroi}(k, :), prd{iroi}(k, :), r2{iroi}(k), exitflg{iroi}(k)] = ...
            dn_fineFit(data, stim1, t1, param, seed, irfType);
    end
end

%% plot fine fit - index
idx = nan(3, 15); midx = zeros(3, 3);

normMax_range = @(x, range) (x - mean(x(range)))./max(x - mean(x(range)));

for iroi = 1 : 3
   for k = 1 : size(roiTS{iroi}, 2)
       
       data = normMax_range(roiTS{iroi}(:, k)', 1 : 200);
       pred = normMax_range(prd{iroi}(k, :), 1 : 200);
%        figure (1), clf
%        plot(data), hold on
%        plot(pred), drawnow, pause(0.5)
 
       idx(iroi, k) = mean(data(700 : 1000) - pred(700 : 1000));
   end
   midx(iroi, 1) = nanmedian(idx(iroi, eccBinIdx{iroi} == 1));
   midx(iroi, 2) = nanmedian(idx(iroi, eccBinIdx{iroi} == 2));
   midx(iroi, 3) = nanmedian(idx(iroi, eccBinIdx{iroi} == 3));
end

figure (1), clf
for k = 1 : 3
    plot(1 : 3, midx(k, :)), hold on
    if k == 3, plot([1, 3], midx(k, [1, 3])), end
end

%ylim([-0.05, 0.2])

%% compute the biphasic fit to V1-V3 electrodes - coarse

irfType = 'biphasic';

for iroi = 1 %: nrois
    for k = 10% 1 : size(roiTS{iroi}, 2)
        data = roiTS{iroi}(200 : 1000, k)';
        [biSeed{iroi}(k, :), biseedR2{iroi}(k)] = dn_gridFit(data, param, stim1(200 : 1000), t1(200 : 1000), irfType);
        k
    end
end

%% compute the biphasic fit to V1 - V3 electrodes - fine fit

% if skip the grid fit ----------------------------
bi = a.dn_biphasic;
%biSeed = bi.biSeed;
irfType = 'biphasic';
roiTS = a.param.roiTS;

for iroi = 1 : 3 %nrois
    iroi 
    for k = 1 : size(roiTS{iroi}, 2)
        k
        %seed = [biSeed{iroi}(k, 1),1, biSeed{iroi}(k, 3 : 4) 0];
        %seed = [biSeed{iroi}(k, :), 0];
        seed = [0.05, 0.8, 0.1, 0.1, 0.001];
        data = roiTS{iroi}(1 : 1000, k)';
        [prm1{iroi}(k, :), prd1{iroi}(k, :), r21{iroi}(k)] = ...
            dn_fineFit(data, stim1(1 : 1000), t1(1 : 1000), param, seed, irfType);
    end
end

%% compute biphasic fit to V1-V3 electrodes - fine fit (with varying "weights" as seed)

weightSeed = linspace(0, 1, 10);

irfType = 'biphasic';
roiTS = a.param.roiTS;

prm1 = {}; prd1 = {}; r21 = {}; 

for k1 = 1 : length(weightSeed)
    seed = [0.05, weightSeed(k1), 0.1, 0.1, 0.001];
    
    for iroi = 1 : 3
        disp([k1, iroi])
        for k = 1 : size(roiTS{iroi}, 2)
            data = roiTS{iroi}(1 : 1000, k)';
            [prm1{iroi}(k1, k, :), prd1{iroi}(k1, k, :), r21{iroi}(k1, k)] = ...
                dn_fineFit(data, stim1(1 : 1000), t1(1 : 1000), param, seed, irfType);
        end  
    end
end


bi_prm = {};
bi_prd = {};
for iroi = 1 : 3
    for k = 1 : size(roiTS{iroi}, 2)
        % find index that corresponds to the best fit:
        idx = find(r21{iroi}(:, k) == max(r21{iroi}(:, k)), 1);
        bi_prm{iroi}(k, :) = squeeze(prm1{iroi}(idx, k, :))';
        bi_prd{iroi}(k, :) = squeeze(prd1{iroi}(idx, k, :))';
    end
end


%% visualize

figure (2), clf

for iroi = 1 : 3 
   for k = 1 : size(roiTS{iroi}, 2)
       data = normMax_range(roiTS{iroi}(1 : 1000, k), 1 : 200);
       pred = normMax_range(bi_prd{iroi}(k, :), 1 : 200);
       
      figure (2), clf
      plot(data), hold on
      plot(pred), drawnow, pause(1)
   end
   
end

%% studying parameters

mweight = [];
sweight = [];

for iroi = 1 : 3
    weight = bi_prm{iroi}(:, 2);
    for ibin = 1 : 3
        mweight(iroi, ibin) = median(weight(eccBinIdx{iroi} == ibin));
        sweight(iroi, ibin, :) = prctile(weight(eccBinIdx{iroi} == ibin), [25, 75]);
    end
end

mmweight = nanmedian(mweight);

figure (3), clf

for k = 1 : 3
   plot(k, mweight(1, k), 'ko', k, mweight(2, k), 'ro', k, mweight(3, k), 'bo', 'markersize', 15), hold on
   plot([k - 0.2, k + 0.2], [mmweight(k), mmweight(k)], 'k-', 'linewidth', 3)
end
set(gca, 'xtick', [1, 2, 3], 'xticklabel', {'low', 'mid.', 'high'}, 'fontsize', 16), xlim([0.5, 3.5]), ylim([0, 1]), box off

%% plot individual time courses

color = {'k', 'b', 'g'};

figure (4), clf

for iroi = 1 : 3
   for k = 1 : size(roiTS{iroi}, 2) 
       subplot_tight(3, 15, (iroi - 1) * 15 + k, 0.02)
       data = normMax_range(roiTS{iroi}(1 : 1000, k), 1 : 200);
       pred = normMax_range(bi_prd{iroi}(k, :), 1 : 200);
       
       plot(data, '-', 'linewidth', 2, 'color', color{eccBinIdx{iroi}(k)}), hold on, plot(pred, 'r', 'linewidth', 2), 
       patch([200, 700, 700, 200], [0, 0, 1, 1], 'k', 'facealpha', 0.1, 'edgecolor', 'w')
       box off, xlim([0, 1000]), ylim([-0.2, 1]), set(gca, 'ytick', [0, 1], 'xtick', '')
       set(gca, 'XAxisLocation', 'origin'), 
       
   end
end
