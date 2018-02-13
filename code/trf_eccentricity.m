% fovea versus periphery

%% load data

dataLoc = fullfile(temporalRootPath, 'results');
fname   = 'params.mat';

a = load(fullfile(dataLoc, fname));
ecog = a.params.ecog;
nroi = length(ecog.roiNm);

roiNm = a.params.ecog.roiNm;

normMax = @(x) x./max(x);

%% plot individual electrode response

labels = ecog.roiElecLabel;
chans  = ecog.roiElecRelAll;
indirsp = ecog.roiElecRsp;

fg = figure;
set(fg, 'position', [1000, 918, 1200, 500]);

for iroi = 1 : nroi
    for k = 1 : length(labels{iroi})
        subplot_tight(nroi, 15, (iroi -1)*15 + k, 0.01)
        plot(indirsp{iroi}(:, k), 'k'), xlim([0, 1200]),ylim([-2, 22]), hold on, axis off
        txt = sprintf('%d', chans{iroi}(k));
        text(700, 12, txt)
    end
end

%% check eccentricity

nroi = 8;

b = load('/Volumes/server/Projects/ECoG/Figures/pRF_data/Power/oneBar/CAR/Subj17_exp_bb.mat');
x0 = b.params.params(:,1);
y0 = b.params.params(:,2);
s0 = b.params.params(:,3)./sqrt(b.params.params(:,5));
for iroi = 1 : nroi
    [x1{iroi}, y1{iroi}, s1{iroi}] = ecogFitPRFPix2Deg(17, x0(chans{iroi}), y0(chans{iroi}), s0(chans{iroi}));
    r1{iroi} = b.params.r(chans{iroi});
    ecc{iroi} = sqrt(abs(x1{iroi}).^2 + abs(y1{iroi}).^2);
end

%% combine later rois:

% Lateral area: the 4th cell, LO
% Ventral area: hV4 and VO, the 5th and 6th
% Dorsal area: V3A and IPS, 7th and 8th cell
t = [1 : length(indirsp{1})]./1000;

cbrsp = {};
cbecc = {};
for k = 1 : 4
    cbrsp{k} = indirsp{k};
    cbecc{k} = ecc{k};
end
cbrsp{5} = [indirsp{5}, indirsp{6}]; cbecc{5} =  [ecc{5}, ecc{6}];
cbrsp{6} = [indirsp{7}, indirsp{8}]; cbecc{6} =  [ecc{7}, ecc{8}];

% average across eccentricities:
low_rsp = []; mid_rsp = []; high_rsp = [];
for k = 1 : 6
    low_rsp(k, :)  = normMax(mean(cbrsp{k}(:, find(cbecc{k} < 5)), 2));
    mid_rsp(k, :)  = normMax(mean(cbrsp{k}(:, find((cbecc{k} >= 5) & cbecc{k} < 10)), 2));
    high_rsp(k, :) = normMax(mean(cbrsp{k}(:, find(cbecc{k} > 10)), 2));
end

%% fit normalization model to the first 700ms of the mean of each eccentricity bin

stim = [zeros(1, 200), ones(1, 500), zeros(1, 504)];
stim1 = [zeros(1, 200), ones(1, 500)];
time = [1 : length(stim1)]./1000;
init = [0.1, 0.1, 3, 0.05, 0.001];
lb   = [0, 0, .5, 0.001, 0];
ub   = [1, 1, 5, 1, 0.5];

m_prm(1, :) = fminsearchbnd(@(x) trf_dcts_finefit(x, normMax(mean(low_rsp(:, 1 : 700))), stim1, time), init, lb, ub);
m_prm(2, :) = fminsearchbnd(@(x) trf_dcts_finefit(x, normMax(nanmean(mid_rsp(:, 1 : 700))), stim1, time), init, lb, ub);
m_prm(3, :) = fminsearchbnd(@(x) trf_dcts_finefit(x, normMax(nanmean(high_rsp(:, 1 : 700))), stim1, time), init, lb, ub);

for k = 1 : 3
    extParam = [m_prm(k, 1), 0, m_prm(k, 2), m_prm(k, 3), m_prm(k, 4), m_prm(k, 5), 1];
    m_prd(k, :) = normMax(trf_dCTSmodel(extParam, stim, t));
end

%trf_dcts_finefit(param, data, stim, t)
%extprm = [param(1), 0, param(2), param(3), param(4), param(5), 1];

%% fit biphasic dCTS to data
fprm = [];
ext_fprm = [];
m_fprm = {};

init = [0.05, 0.9, 0.07, 0.1, 0.001];
lb   = [0, 0, 0, 0.001, 0];
ub   = [1, 1, 1, 1, 0.5];

for k = 1 : 6
   for k1 = 1 : size(cbrsp{k}, 2)
       data = normMax(cbrsp{k}(:, k1));
       fprm(k, k1, :) = fminsearchbnd(@(x) trf_dcts_biphasicFit(x, data', stim, t), init, lb, ub);
       ext_fprm(k, k1, :) = [fprm(k, k1, 1),fprm(k, k1, 2), fprm(k, k1, 3), 2, fprm(k, k1, 4), fprm(k, k1, 5), 1];
       m_fprd{k}(k1, :) = normMax(trf_dCTSmodel(squeeze(ext_fprm(k, k1, :)), stim, t));
   end
end

%% plot biphaisc fit

figure (5), clf
for k = 1 : 6
    for k1 = 1 :  size(cbrsp{k}, 2)
        subplot_tight(6, 15, (k - 1) * 15 + k1, 0.02),
        plot(normMax(cbrsp{k}(:, k1)), 'k-'), hold on
        plot(abs(m_fprd{k}(k1, :)), 'r-')
        axis off, axis tight
    end
end

%% plot fitted parameters

lm = nan(1, 6);
mm = nan(1, 6);
hm = nan(1, 6);
figure (6), clf

for k = 1 : 6
    subplot(6, 1, k)
    
    lm(k) = mean(fprm(k, find(cbecc{k}<=5), 2));
    hm(k) = mean(fprm(k, find(cbecc{k} > 10), 2));
    
    plot(1, mean(fprm(k, find(cbecc{k}<=5), 2)), 'ro'), hold on
    if ~ismember(k, [3, 4])
        mm(k) = mean(fprm(k, find((cbecc{k}>5 & cbecc{k} <= 10)), 2));
        plot(2, mean(fprm(k, find((cbecc{k}>5 & cbecc{k} <= 10)), 2)), 'bo'),
    end
     plot(3, mean(fprm(k, find(cbecc{k}>10), 2)), 'ko'),
     xlim([0.5, 3.5])
end

figure (7), clf
plot([0.8, 1.2], [mean(lm(1 : 3)), mean(lm(1 : 3))], 'r-'); hold on
plot([1.8, 2.2], [nanmean(mm(1 : 3)), nanmean(mm(1 : 3))], 'r-'), 
plot([2.8, 3.2], [mean(hm(1 : 3)), mean(hm(1 : 3))], 'r-')

p1 = plot(1, lm(1 : 3), 'ko', 2, mm(1 : 3), 'ko', 3, hm(1 : 3), 'ko');

set(p1, 'markersize', 5, 'markerfacecolor', 'k', 'markeredgecolor', 'r')
xlim([0.5, 3.5]), box off, ylim([0, 1]), set(gca, 'xtick', 1 : 3)


%% plot averaged eccentricity response in each roi

figure (2), clf
for k = 1 : 3
    subplot(1, 3, k)
    patch([0.7, 1, 1, 0.7], [0, 0, 1, 1], [0.9, 0.9, 0.9]), hold on
    patch([0.2, 0.7, 0.7, 0.2], [0, 0, 1, 1], [1, 1, 1])
    ax = gca; ax.XAxisLocation = 'origin';
    set(gca, 'xtick', [0.2, 0.7, 1]), set(gca, 'ytick', [0, 1])
    %plot(t, m_prd(k, :), 'r:')
end
subplot(1, 3, 1), plot(t, normMax(nanmean(low_rsp)), 'k-'), axis tight, box off, hold on,
subplot(1, 3, 2), plot(t, normMax(nanmean(mid_rsp)), 'k-'), axis tight, box off, hold on
subplot(1, 3, 3), plot(t, normMax(nanmean(high_rsp)), 'k-'), axis tight, box off, hold on

%% fit the first 700ms of each individual data

init = [0.1, 0.1, 2, 0.05, 0.001];

for k = 1 : 6
    if k < 4, init = [0.07, 0.07, 1.5, 0.1, 0.001];
    else, init = [0.2, 0.2, 2.5, 0.01, 0.001]; end
   for k1 = 1 : size(cbrsp{k}, 2)
       prm{k}(k1, :) = fminsearchbnd(@(x) trf_dcts_finefit(x, normMax(cbrsp{k}(1 : 700, k1))', stim1, time), init, lb, ub);
       extParam = [prm{k}(k1, 1), 0, prm{k}(k1, 2), prm{k}(k1, 3), prm{k}(k1, 4), prm{k}(k1, 5), 1];
       prd{k}(k1, :) = normMax(trf_dCTSmodel(extParam, stim, t));
   end
end

%% compute the index

idx = nan(6, 15);

for k = 1 : 6
   for k1 = 1 : size(cbrsp{k}, 2)
      idx(k, k1) = mean(cbrsp{k}(700 : 1204, k1)'- prd{k}(k1, 700 : 1204));
   end
   lowidx(k) = mean(idx(k, find(cbecc{k} < 5)));
   mididx(k) = mean(idx(k, find((cbecc{k} >= 5) & cbecc{k} < 10)));
   highidx(k)= mean(idx(k, find(cbecc{k} >= 10)));
end

figure (4), clf
for k = 1 : 6
   subplot(1, 6, k)
   bar(1, lowidx(k)), hold on
   bar(2, mididx(k))
   bar(3, highidx(k)), xlim([0.5, 3.5]), ylim([-0.1, 1]), box off
end
% plot(lowidx, 'ko-'), hold on
% plot(mididx, 'bo-')
% plot(highidx, 'ro-')



%% visualize the fit
figure (3), clf
for k = 1 : 6
    for k1 = 1 :  size(cbrsp{k}, 2)
        subplot_tight(6, 15, (k - 1) * 15 + k1, 0.02),
        plot(normMax(cbrsp{k}(:, k1)), 'k-'), hold on
        plot(prd{k}(k1, :), 'r-')
        axis off, axis tight
    end
end

%% fit the entire time series of each individual data




