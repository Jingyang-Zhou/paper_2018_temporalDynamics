% dn_mkFigure3

%% load data

dataLoc = fullfile(temporalRootPath, 'data');
fname = 'dn_data.mat';

a    = load(fullfile(dataLoc, fname));
dn = a.dn;

nboots = 100;
nrois  = 6;

param = dn.param;
%param = param(dn.exitflg, :);

rois = {'V1', 'V2', 'V3', 'LA', 'VA', 'DA'};

t = [1 : size(a.bsData, 3)]/1000;

%% summarize parameters

mparam = [];
sparam = [];

for k = 1 : nrois
    idx = (k - 1) * nboots + 1 : k * nboots;
    mparam(k, :) = mean(param(idx, 1 : 4));
    %sparam(k, :, :) = prctile(param(idx, 1 : 4), [25, 75], 1);
    sparam(k, :, :) = std(param(idx, :))./2;
end

figure (1), clf
for k = 1 : 4
    subplot(1, 4, k)
    plot(1 : nrois, mparam(:, k), 'ko'), hold on
   for iroi = 1 : nrois
      % plot([iroi, iroi], squeeze(sparam(iroi, :, k)), 'k-')
      plot([iroi, iroi], [mparam(iroi, k) - sparam(iroi, k),mparam(iroi, k) + sparam(iroi, k)], 'k-')
   end
   xlim([0.5, nrois + 0.5])
   set(gca, 'xtick', [1 : 6], 'xticklabel', rois), box off
end

subplot(1, 4, 3), ylim([0, 4])
subplot(1, 4, 4), ylim([0, 0.2])

%% plot summary metrics

for k = 1 : nrois
    idx = (k - 1) * nboots + 1 : k * nboots;
    mt2pk(k) = median(dn.derivedPrm.t2pk(idx));
    st2pk(k, :) = prctile(dn.derivedPrm.t2pk(idx), [25, 75]);
    
    masymp(k) = median(dn.derivedPrm.r_asymp(idx));
    sasymp(k, :) = prctile(dn.derivedPrm.r_asymp(idx), [25, 75]);
end


figure (2), clf
subplot(1, 2, 1)
plot(mt2pk, 'ko', 'markersize', 10, 'markerfacecolor', 'k'), hold on
for k = 1 : nrois
   plot([k, k], [st2pk(k, :)], 'k-') 
end
xlim([0.5, nrois + 0.5]), ylim([0, 200]), box off

subplot(1, 2, 2)
plot(masymp, 'ko', 'markersize', 10, 'markerfacecolor', 'k'), hold on
for k = 1 : nrois
   plot([k, k], [sasymp(k, :)], 'k-') 
end
xlim([0.5, nrois + 0.5]), ylim([0, 0.25]), box off

%% plot ECoG time courses and model fit


% plot data
figure (3), clf

for k = 1 : length(rois)
    subplot(2, 3, k), 
    m = squeeze(median(a.bsData(k, :, :), 2));
    s = squeeze(prctile(a.bsData(k, :, :), [25, 75], 2));
    
    shadedErrorBar(t, m, [m'-s(1, :); s(2, :) - m'], 'k-'), hold on, 
    patch([0.2, 0.7, 0.7, 0.2], [0, 0, 17, 17], 'k', 'facealpha', 0.05)
    
    % plot model prediction
    idx = (k - 1) * 100 + 1 : k * 100;
    m1 = median(dn.prd(idx, :).*dn.param(idx, 6));
    %s1 = prctile(dn.prd(idx, :).*dn.param(idx, 6), [25, 75]);
    %shadedErrorBar(t, m1, [m1-s1(1, :); s1(2, :) - m1], 'r-')
    plot(t, m1, 'r-', 'linewidth', 2), 
    set(gca, 'xaxislocation', 'origin', 'xtick', [0.2, 0.7], 'xticklabel', [0, 0.5])
    
    xlim([0, max(t)]), ylim([-1, 17]), box off
end
