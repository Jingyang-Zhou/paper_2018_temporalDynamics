% trf_contrast_modulation

dataLoc = fullfile(dn_ECoG_RootPath, 'data');
fName   = 'Figure1_ACE_Data.xlsx';
a = xlsread(fullfile(dataLoc, fName));

cell1 = a(2 : 27, 2 : 11);

%% extract numbers

% make stimulus
stim = [ones(1, 200), zeros(1, 300)];

%contrast_levels = a(1, 2 : end);
contrast_levels = 0 : 0.1 : 0.9;
cellRsp{1}  = a(2 : 27, 2 : 11); 
time{1}     = a(2 : 27, 1);

cellRsp{2}  = a(30 : 50, 2 : 11);
time{2}     = a(30 : 50, 1);

cellRsp{3}  = a(53 : 73, 2 : 11);
time{3}     = a(53 : 73, 1);

% normalize the response to the max
normMax = @(x) x./max(x(:));

for k = 1 : 3
   cellRsp{k} = normMax(cellRsp{k});
end

% visualize the response to different contrast levels
figure (1), clf
for k = 1 : 3
    subplot(1, 3, k)
    set(gca, 'ColorOrder', parula(15)); hold on
    plot(time{k}, cellRsp{k}), xlabel('time(ms)'), ylabel('response(arbitrary unit)'), axis tight
end

%% fit model

init = [0.05, 0.05, 2, 0.01, 0.03];
low  = [0, 0, 0, 0, 0];
high = [1, 1, 10, 1, 0.1];
prms = [];

pred = []; 
for k = 1 : 3
    prms(k, :)   = fminsearchbnd(@(x) dn_fitModel_contrast(x, cellRsp{k}', time{k}, contrast_levels, stim), init, low, high);
    [~, pred{k}] = dn_fitModel_contrast(prms(k, :), cellRsp{k}, time{k}, contrast_levels, stim);
    tmp = pred{k};
    pred{k} = pred{k} / max(tmp(:));
end

%% COMPUTE VARIANCE EXPLAINED FOR EACH CELL

for k = 1 : 3
    tmp = pred{k}(:, time{k});
    model = tmp(:);
    tmp = cellRsp{k}';
    data  = tmp(:);
    top = sum((model - data).^2);
    bottom = sum(data.^2);
    r2(k) = 1-top /bottom
end

%% plot data versus model fit

figure (2), clf
for k = 1 : 3
    subplot(2, 3, k)
    set(gca, 'ColorOrder',copper(9)); hold on
    plot(time{k}, cellRsp{k}(:, 2 : end), '-', 'linewidth', 2), axis tight, axis square, box off, 
    set(gca, 'ytick', [0, 1])
   
    subplot(2, 3, 3 + k)
    set(gca, 'ColorOrder', copper(9)); hold on
    plot(time{k}, pred{k}(2 : end, time{k}), '-'), axis tight, axis square, box off, set(gca, 'ytick', [0, 1])
end

 subplot(2, 3, 1)
 col = bone(12);

figure (3), clf
for k = 1 : size(prms, 2) - 1
    subplot(2, 2, k)
    set(gca, 'ColorOrder', bone(12)); hold on
    plot([1, 1, 1], prms(:, k), 'o')
end

%% plot maximum heights and value of the maximum heights

d_max_time = {};

for k = 1 : 3
   for k1 = 1 : size(cellRsp{1}, 2)
       d_max_time{k}(k1) = find(cellRsp{k}(:, k1) == max(cellRsp{k}(:, k1)), 1);
       d_max_height{k}(k1) = cellRsp{k}(d_max_time{k}(k1), k1);
       
       m_max_time{k}(k1) = find(pred{k}(k1, :) == max(pred{k}(k1, :)), 1);
       m_max_height{k}(k1) = pred{k}(k1, m_max_time{k}(k1));
   end
   d_max_time{k} = time{k}(d_max_time{k});
   m_max_time{k} = time{k}(m_max_time{k});
end
figure (4), clf
col = copper(9);
lim = {[65, 120], [45, 83], [55, 115]};
for k = 1 : 3
    for k1 = 2 : size(m_max_time{1}, 1)
        subplot(2, 3, k), plot(d_max_time{k}(k1), m_max_time{k}(k1), 'o', 'markerFaceColor', col(k1-1, :), 'markeredgecolor', 'k', 'markersize', 10),
        hold on,
        subplot(2, 3, k + 3), plot(d_max_height{k}(k1), m_max_height{k}(k1), 'o', 'markeredgecolor', 'w', 'markerfacecolor', col(k1-1, :), 'markersize', 12), hold on,
    end
    subplot(2, 3, k),  plot([0, 120], [0, 120], 'k:'),  xlim(lim{k}), ylim(lim{k}), axis square
    set(gca, 'xtick', lim{k}(1) : 30 : lim{k}(2), 'ytick', lim{k}(1) : 30 : lim{k}(2)),
    subplot(2, 3, k + 3),  plot([0, 1], [0, 1], 'k:'), axis square
end

figure (5), clf
set(gca, 'ColorOrder', copper(9)); hold on
plot(time{1}, cellRsp{1}(:, 2: 10)), axis square
for k = 9 : -1 : 1
   txt{k} = sprintf('ctrst. = %d', k *10); 
end
legend(txt)

%% load ECoG data

dataLoc = fullfile(temporalRootPath, 'data');
fName   = 'dn_data.mat';
a       = load(fullfile(dataLoc, fName));
dn      = a.dn;

ecog.msummaryPrm =[median(dn.derivedPrm.t2pk(1 : 100)), median(dn.derivedPrm.r_asymp(1 : 100))];
ecog.msummaryPrm =[prctile(dn.derivedPrm.t2pk(1 : 100), [25, 75]); prctile(dn.derivedPrm.r_asymp(1 : 100), [25, 75])];

%% compute summary metrics and plot it against the ECoG data

summary = dn_computeDerivedParams(prms, 'uniphasic');

figure (6), clf
yyaxis left, plot(1, summary.t2pk, 'ko'), hold on, plot(1, ecog.msummaryPrm(1), 'ro'), ylim([0, 250])
yyaxis right, plot(2, summary.r_asymp, 'ko'), plot(2, ecog.msummaryPrm(2), 'ro'), ylim([0, 0.5])

xlim([0.5, 2.5])

%%

%% MAKE A STIMULUS

t    = 0.001 : 0.001 : .2;
stim = zeros(size(t));
stim(1 : 200) = 1;

% make contrast stimulus

nContrast = 10;
contrast  = linspace(0.1, 1, nContrast);
cStim     = zeros(nContrast, length(t));
for k = 1 : nContrast
   cStim(k, :) = stim .* contrast(k);
end

%% MAKE DN MODEL PREDICTIONS

prd = [];

param = [0.02, 0, 0.03, 2, 0.2, 0, 1];

for k = 1 : nContrast
   prd(k, :) = dn_DNmodel(param, cStim(k, :), t);
end

figure (1), clf
subplot(121), set (gca, 'colororder', copper(nContrast)), hold on, plot(t, cStim, 'linewidth', 2)
subplot(122), set (gca, 'colororder', copper(nContrast)), hold on, plot(t, prd, 'linewidth', 2), 
xlabel('time (s)'), set(gca, 'ytick', '', 'fontsize', 14),

%% GENERATE A NAKA-RUSHTON EQUATION

prm = [];

prm(1, :) = [0.0395, 0, 0.0397, 2.2829, 0.0266, 0, 1];
prm(2, :) = [0.0436, 0, 0.1049, 1.6136, 0.0027, 0, 1];
prm(3, :) = [0.0337, 0, 0.1373, 3.5056, 0.016,  0, 1];

cprd = [];

for k = 1 : nContrast
   for icell = 1 : 3
      cprd(icell, k, :) = dn_DNmodel(prm(icell, :), cStim(k, :), t); 
   end
end

figure (2), clf
for k = 1 : 3
   subplot(2, 3, k), set(gca, 'colororder', copper(nContrast)), hold on
   cprd(k, :, :) = cprd(k, :, :) ./max(max(cprd(k, :)));
   plot(t, squeeze(cprd(k, :, :)), 'linewidth', 2)
end

%% FIT THE SHIFTED AND SCALED COPY MODEL

init = [1, 0.04];
lb   = [0, 0];
ub   = [2, 0.15];

% compute average shape
for k = 1 : 3
    tmp = squeeze(cprd(k,  :, :));
    for k1 = 1 : nContrast
       tmp(k1, :) = tmp(k1, :) ./ max(tmp(k1, :));
    end
   aveshape(k, :) =  mean(tmp);
end

for k = 1 : 3
    example = squeeze(cprd(k, 8, :));
    %example = aveshape(k, :)';
   for k1 = 1 : nContrast
       sprm(k, k1, :) = fminsearchbnd(@(x) fitShiftedScaledCopy(example, x, squeeze(cprd(k, k1, :))), init, lb, ub);
       sprd(k, k1, :) = computeShiftedScaledCopy(example, squeeze(sprm(k, k1, :))); 
   end
   % compute R2
   tmp1 = sprd(k, :);
   tmp2 = cprd(k, :);
   top = sum((tmp1 - tmp2).^2);
   bottom = sum(tmp2.^2);
   r2(k) = 1-top/bottom
end

figure (2), 
for k = 1 : 3
   subplot(2, 3, k + 3), cla,  set(gca, 'colororder', copper(nContrast)), hold on
   plot(t, squeeze(sprd(k, :, :))', 'linewidth', 2), box off, ylim([0, 1])
   
end


