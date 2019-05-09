% dn_fitDN2ContrastSUA
%
% DESCRIPTION -------------------------------------------------------------
% In this file we fit the DN model to single cells' response to a static
% stimulus of different contrast levels. The DN model is fit to all
% response time courses to 10 different contrast levels concurrently.

%% SAVE FIGURE KNOB

saveFigure = 0;

%% LOAD ALBRECHT AND GEISLER DATA

dataLoc = fullfile(dn_ECoG_RootPath, 'data');
fName   = 'Figure1_ACE_Data.xlsx';
a       = xlsread(fullfile(dataLoc, fName)); % concatenated time course (over 3 cells) x contrast levels

%% VISUALIZE THE ORIGINAL DATA a
% 
% figure (100), clf, set(gca, 'colororder', copper(11)), hold on
% plot(a, 'linewidth', 2), axis tight, set(gca, 'fontsize', 14)
% xlabel('first dimension'), ylabel('second dimension')

%% EXTRACT AND DERIVE PARAMETERS

% MAKE STIMULUS -----------------------------------------------------------
stim = [ones(1, 200), zeros(1, 300)];
t    = 1 : 1 : length(stim);

% EXTRACT RESPONSE FOR EACH CELL ------------------------------------------
contrast_levels = a(1, 2 : end);
cellRsp{1}  = a(2 : 27, 2 : 11)./100; % the time courses have maximal repsonse 100, and we normalize to maximal response 1 here
time{1}     = a(2 : 27, 1); % in unit of millisecond                                                                                                                                           []p

cellRsp{2}  = a(30 : 50, 2 : 11)./100; % response from low to hight
time{2}     = a(30 : 50, 1);

cellRsp{3}  = a(53 : 73, 2 : 11)./100;
time{3}     = a(53 : 73, 1);

%% FIT DN MODEL TO THE CONTRAST DATA

init = {[0.1, 0.1, 2.5, 0.1, 0.04], [0.1, 0.1, 4, 0.01, 0.035], [0.1, 0.1, 3, 0.1, 0.04]};
low  = [0, 0, 0, 0, 0];
high = [0.8, 1, 10, 1, 0.1];
sua  = [];

pred = [];
for iCell = 1 : 3
    sua.ctrprms(iCell, :) = fminsearchbnd(@(x) dn_fitModel_contrast(x, cellRsp{iCell}', time{iCell}, contrast_levels, stim), init{iCell}, low, high);
    [~, sua.ctrprd{iCell}] = dn_fitModel_contrast(sua.ctrprms(iCell, :), cellRsp{iCell}', time{iCell}, contrast_levels, stim);
    sua.ctrprd{iCell} = sua.ctrprd{iCell}./max(sua.ctrprd{iCell}(:));
end

%% COMPUTE TIME TO PEAK AND PEAK RESPONSE LEVEL

dt_t2pk = []; prd_t2pk = []; dt_pklevel = []; prd_pklevel = [];

for icell = 1 : 3
    % DATA
    [dt_pklevel(icell, :), dt_t2pk(icell, :)] = max(cellRsp{icell});
    dt_t2pk(icell, :) = time{icell}(dt_t2pk(icell, :));
    % PREDICTION
    [prd_pklevel(icell, :), prd_t2pk(icell, :)] = max(sua.ctrprd{icell}(:, 1 : 200)');
    prd_t2pk(icell, :) = t(prd_t2pk(icell, :));
end

%% VISUALIZE SUR (SINGLE UNIT RESPONSE) TO DIFFERENT CONTRAST LEVELS

contrast = 0 : 10 : 90;

fig_start = [1, 3, 5]; title_txt = {'CELL 1', 'CELL 2', 'CELL 3'};

fg1 = figure (1); clf
for iCell = 1 : 3
    % PLOT DATA -----------------------------------------------------------
    subplot(5, 6, fig_start(iCell) : fig_start(iCell) + 1), set(gca, 'ColorOrder', copper(10)); hold on
    t0 = min(time{iCell}); t1 = max(time{iCell});
    % patch([t0, t1, t1, t0], [0, 0, 1, 1], 0.9 * ones(1, 3))
    plot(time{iCell}, cellRsp{iCell}, 'linewidth', 2), xlabel('time(ms)'),
    axis tight, set(gca, 'fontsize', 14, 'ytick', [0, 1]), title(title_txt{iCell})
    if iCell == 1, ylabel('response(arbitrary unit)'), end, 
    
    % PLOT MODEL PREDICTION -----------------------------------------------
    subplot(5, 6, 6 + fig_start(iCell) : fig_start(iCell) + 7), set(gca, 'ColorOrder', copper(10)); hold on
    plot(t, sua.ctrprd{iCell}, 'linewidth', 2), xlabel('time (ms)'), xlim([t0, t1])
    set(gca, 'fontsize', 14, 'ytick', [0, 1]), title(title_txt{iCell}), 
    
    subplot(5, 6, 12 + fig_start(iCell)),
    plot(contrast(2 : end), dt_t2pk(iCell, 2 : end), 'ko', 'markerfacecolor', 'k', 'markersize', 7), hold on
    plot(contrast(2 : end), prd_t2pk(iCell, 2 : end), 'r-', 'linewidth', 3)
    axis square, box off, axis tight, xlim([0, 100]), 
    xlabel('% contrast'), set(gca, 'ytick', [dt_t2pk(iCell, end), dt_t2pk(iCell, 2)]), set(gca, 'xtick', [10, 90], 'fontsize', 14)
    
    subplot(5, 6, 13 + fig_start(iCell))
    plot(contrast(2 : end), dt_pklevel(iCell, 2 : end), 'ko', 'markerfacecolor', 'k', 'markersize', 7), hold on
    plot(contrast(2 : end), prd_pklevel(iCell, 2 : end), 'b-', 'linewidth', 3)
    axis square, box off, axis tight, xlim([0, 100]),   set(gca, 'xtick', [10, 90], 'ytick', [0, 1], 'fontsize', 14), 
end

fg1.Position = [0, 1000, 1200, 1200];


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

subplot(5, 6, [19 : 20, 25 : 26]), set (gca, 'colororder', copper(nContrast)), hold on, plot(t, prd, 'linewidth', 2), 
xlabel('time (s)'), set(gca, 'ytick', '', 'fontsize', 14), title('Non-converging rsp. at stimulus offset')

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

for k = 1 : 3
   subplot(5, 6, 21+k), set(gca, 'colororder', copper(nContrast)), hold on
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
   r2(k) = 1-top/bottom;
end

for k = 1 : 3
   subplot(5, 6, 27 + k), cla,  set(gca, 'colororder', copper(nContrast)), hold on
   plot(t, squeeze(sprd(k, :, :))', 'linewidth', 2), box off, ylim([0, 1])
end





