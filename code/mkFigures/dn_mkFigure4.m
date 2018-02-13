% dn_mkFigure 4

%% load data

dataLoc = fullfile(temporalRootPath, 'data');
fName   = 'dn_data.mat';

a       = load(fullfile(dataLoc, fName));
param   = a.param;

roiNm = param.roiNm;
nrois = length(roiNm);
roiTS = param.roiTS;

bi = a.dn_biphasic;

normMax = @(x) x./max(x);

%% plot panel A

% average data according to eccentricities bins, and plot the difference
biidx = a.param.eccBinIdx;
eccBinTS = {};
t = [1 : length(roiTS{1})]./1000;
stim = zeros(1, length(t)); stim([201 : 700]) = 1;

for iroi = 1 : nrois
    for k = 1 : 3
        idx = find(biidx{iroi} == k);
        if iroi ~= 1
            eccBinTS{k} = [eccBinTS{k}, normMax(roiTS{iroi}(:, idx))];
        else
            eccBinTS{k} = normMax(roiTS{iroi}(:, idx));
        end
    end
end

% % fit the monophasic model to the first 700ms of the stimulus
% t1 = 0.001 : 0.001 : 0.7;
% stim1 = zeros(1, length(t1)); stim1(201 : end) = 1;
% % seed = [0.1, 0.9, 3, 0.001, 0];
% for k = 1 : 3
%     meccBinTS(k, :) = normMax(mean(eccBinTS{k}, 2));
%     % grid fit
%     [modelSeed, seedR2] = dn_gridFit(meccBinTS(:, 1 : 700), param, stim1, t1, 'uniphasic');
%     % fine fit
%     [prm(k, :), prd0(k, :), ~] = dn_fineFit(meccBinTS(k, 1 : 700), stim1, t1, param, [modelSeed(k, :), 0], 'uniphasic');
%     prd(k, :) = normMax(dn_DNmodel([prm(k, 1), 0, prm(k, 2 : 5), 1], stim, t));
% end


% figure (1), clf
% t_postStim = 700 : 1000;
% for k = 1 : 3
%     %meccBinTS = normMax(mean(eccBinTS{k}, 2));
%     max_postStim = find(meccBinTS(k, t_postStim) == max(meccBinTS(k, t_postStim))) + 700;
%     subplot(1, 3, k)
%     patch([0.2, 0.7, 0.7, 0.2], [0, 0, 1, 1],  1*ones(1, 3)), hold on
%     patch([0.7, 1, 1, 0.7], [0, 0, meccBinTS(k,max_postStim), meccBinTS(k,max_postStim)], 'y', 'facealpha', 0.3)
%     plot(t, meccBinTS(k, :), 'k-', 'linewidth', 2), ylim([-0.1, 1]), xlim([0, 1.204]), box off
%     %plot(t, prd(k, :), 'r:')
%     set(gca, 'ytick', [0, 1]), set(gca, 'xtick', [0.2, 0.7], 'xticklabel', [0, 0.5])
%     set(gca, 'XAxisLocation', 'origin')
%
% end

%% plot panel B

%% plot individual fit to the first 700ms of the stimulus

ts_color = ['k', 'b', 'g'];

figure(2), clf
for iroi = 1 : nrois
    for k = 1 : size(roiTS{iroi}, 2)
        subplot_tight(6, 15, (iroi - 1) * 15 + k, 0.01)
        plot(t, normMax(roiTS{iroi}(:, k)), '-', 'color', ts_color(biidx{iroi}(k)), 'linewidth', 3), hold on,
        plot(t, normMax(bi.ORIprd{iroi}(k, :)), 'r-', 'linewidth', 2)
        % outline the region that we did our analysis on
        patch([0.7, 1, 1, 0.7], [-0.2, -0.2, 1, 1], 'r', 'facealpha', 0.2)
        set(gca, 'XAxisLocation', 'origin')
        ylim([-0.2, 1]), xlim([0, 1.204]), box off, set(gca, 'xtick', [0.2, 0.7, 1], 'xticklabel', '')
        set(gca, 'ytick', [0, 1]), if k ~= 1, set(gca, 'yticklabel', ''), end
    end
end
%% plot panel C

% plot the offset response index
for iroi = 1 : nrois
    % compute the mean difference between prediction and measurement after
    % stimulus offset
    
    for k = 1 : size(roiTS{iroi}, 2)
        rs = normMax(roiTS{iroi}(:, k));
        pd = normMax(bi.ORIprd{iroi}(k, :)');
        ORIidx{iroi}(k) = mean(rs(700 : 1000) - pd(700 : 1000));
    end
    %     ts = normMax(roiTS{iroi}())
    %
    %     diff = normMax(roiTS{iroi}(700 : 1000, :))' - normMax(bi.ORIprd{iroi}(:, 700 : 1000));
    %     ORIidx{iroi} = max(diff, [], 2);
end

% average across different eccentricity bins
for iroi = 1 : nrois
    for k = 1 : 3
        idx = find(biidx{iroi} == k);
        binned_ROIidx{iroi, k} =  ORIidx{iroi}(idx);
        mbinned_ROIidx(iroi, k) = mean(binned_ROIidx{iroi, k})
        sbinned_ROIidx(iroi, k) = std(binned_ROIidx{iroi, k});
    end
end

figure (3), clf

eccColor = ['k', 'b', 'g']

for iroi = [1, 2, 5, 6]
    subplot(1, 6, iroi)
    m = mbinned_ROIidx(iroi, :);
    s = sbinned_ROIidx(iroi, :);
    plot(m, 'ko:', 'markerfacecolor', 'k', 'markersize', 10), hold on,
    for k = 1 : 3,
        plot([k, k], [m(k) - s(k), m(k) + s(k)], 'k-'),
    end
    box off, xlim([0.5, 3.5])
end

for iroi = [3, 4]
    subplot(1, 6, iroi)
    m = [mbinned_ROIidx(iroi, 1),0, mbinned_ROIidx(iroi, 3)];
    s = [sbinned_ROIidx(iroi, 1),0, sbinned_ROIidx(iroi, 3)];
    plot([1, 3], m([1, 3]), 'ko:', 'markerfacecolor', 'k', 'markersize', 10), hold on
    for k = [1, 3]
        plot([k, k], [m(k) - s(k), m(k) + s(k)], 'k-')
    end
    box off, xlim([0.5, 3.5])
end

for k = 1 : 6
    subplot(1, 6, k)
    set(gca, 'xtick', [1 : 3], 'xticklabel', {'low', 'mid', 'high'}), ylim([-0.07, 0.25])
    title(roiNm{k}), ylabel('index')
end


%% plot panel D : visualize time courses and biphasic model fit

figure (4), clf
for iroi = 1 : nrois
    for k = 1 : size(roiTS{iroi}, 2)
       subplot_tight(6, 15, (iroi - 1) * 15 + k, 0.01)
        plot(normMax(roiTS{iroi}(:, k))), hold on
        plot(normMax(bi.biprd{iroi}(k, :)), 'r'), box off, axis tight, ylim([-0.2, 1])
    end
end

%% plot panel E

prm = bi.biprm;

for iroi = 1 : nrois
    for k = 1 : 3
        idx = find(param.eccBinIdx{iroi} == k);
        mprm(iroi, k) = median(prm{iroi}(idx, 2));
        sprm(iroi, k) = std(prm{iroi}(idx, 2));
    end
end

figure

plot(mprm(2, :), 'o:'), xlim([0.5, 3.5])

%% visualize model fit

bi = a.dn_biphasic;
normMax = @(x) x./max(x);

figure (1), clf
for k = 1 : 12
    subplot(12, 1, k)
    plot(normMax(param.roiTS{1}(:, k))), hold on
    plot(normMax(bi.biprd{1}(k, :)), 'r')
end
