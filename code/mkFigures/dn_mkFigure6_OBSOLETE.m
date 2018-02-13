% dn_mkFigure6

%% useful functions

normMax_range = @(x, range) (x - mean(x(range)))./max(x - mean(x(range)));
normMax       = @(x) x./max(x);

plot_gap = 0.03;

%% load data:

% load single unit data ---------------------------------------------------
dataLoc1 = fullfile(temporalRootPath, 'data', 'ECoG'); fName1 = 'figure4Data.xlsx'; fName2 = 'Figure1_ACE_Data.xlsx';
% load data 
a      = xlsread(fullfile(dataLoc1, fName1));
% extract data
t_sg   = a(:, 1)./1000;
dt_sg  = a(:, 2 : end);
nCells = size(dt_sg, 2);

% load the second set of data
% b      = xlsread(fullfile(dataLoc1, fName2));
% dt_sg1{1} = b(2 : 27, 11);  t_sg1{1} = b(2 : 27, 1)./1000;
% dt_sg1{2} = b(30 : 50, 11); t_sg1{2} = b(30 : 50, 1)./1000;
% dt_sg1{3} = b(53 : 73, 11); t_sg1{3} = b(53 : 73, 1)./1000;

% load multi-unit data ----------------------------------------------------
fName3 = 'CONTEXT.mat';
% load data
c = load(fullfile(dataLoc1, fName3));
% extract data
elecIdx = c.DETS(:, 1); sessIdx = c.DETS(:, 2); % electrode and session index

dt_mua{1} = c.MUA(elecIdx == 6, 2 : end); dt_mua{2} = c.MUA(elecIdx == 7, 2 : end);
dt_lfp{1} = c.LFP(elecIdx == 6, 2 : end); dt_lfp{2} = c.LFP(elecIdx == 7, 2 : end);

srate = c.FsD;

% load ECoG data ----------------------------------------------------------
fName4 = 'dn_data.mat';
% load data
d = load(fullfile(fileparts(dataLoc1), fName4));
dt_ecog = d.bsData;

%% PART 1: FIT DN MODEL TO THE SINGLE UNIT DATA

% fit model to response from a couple cells, and estimate time constants

init    = [0.04, 0.05, 1.8, 0.03, 0.01]; % 'tau1',  'tau2', 'n', 'sigma', 'shift'
lb      = [0.01,0.01,0,0,0];
ub      = [1, 1, 10, 1, 1];

cells2Analyze = [10, 11, 12];
stim_sg = ones(1, length(t_sg)); % stimulus on

sg   = struct(); extended_stim = padarray(stim_sg, [0, 600], 'post'); extended_t_sg = [1 : length(extended_stim)]./1000;

for k = cells2Analyze
    sg.param(k - 9, :) = fminsearchbnd(@(x) dn_computeFineFit(x, dt_sg(:, k)', stim_sg, t_sg, 'uniphasic'), init, lb, ub);
    sg.prd(k - 9, :)   = normMax_range(dn_DNmodel([sg.param(k - 9, 1), 0, sg.param(k - 9, 2 : end), 1], extended_stim, t_sg), [1 : 20]);
end

% summary metrics
sg.metrics = dn_computeDerivedParams(sg.param, 'uniphasic');

% sg1 = [];
% 
% % prd_time = 
% 
% for k = 1 : 3
%     stim_sg1{k} = ones(1, length(t_sg1{k}));
%     sg1.param(k, :) = fminsearchbnd(@(x) dn_computeFineFit(x, dt_sg1{k}', stim_sg1{k}, t_sg1{k}, 'uniphasic'), init, lb, ub);
%     %sg1.prd(k, :)   = normMax(dn_DNmodel([sg1.param(k, 1), 0, sg1.param(k, 2 : end), 1], stim_sg1{k}, t_sg1{k}));
% end
% 
% sg.metrics = dn_computeDerivedParams(sg1.param, 'uniphasic');

%% PART 1: VISUALIZATION

figure (1), clf
for k = 1 : 3
   subplot_tight(4, 4, k, plot_gap), cla
   stem(t_sg, dt_sg(:, k+9), 'k-', 'markersize', 1), hold on, 
   plot(extended_t_sg, sg.prd(k, :), 'r-', 'linewidth', 2)
   patch([0, 0.2, 0.2, 0], [0, 0, 1, 1], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'w'), 
   xlim([-0.2, 0.8]), box off, ylim([-0.2, 1]), axis square, ax = gca; ax.XAxisLocation = 'origin';
   set(gca, 'xtick', [0, 0.2], 'ytick', [0, 1], 'fontsize', 14), title('PSTH (single unit)')
end

subplot_tight(4, 4, 4, plot_gap), cla, grayOrder = [0.9, 0.6, 0.1];
for k = 1 : 3
    yyaxis left
    plot(1, sg.metrics.t2pk(k), 'ko', 'markerfacecolor', grayOrder(k)*ones(1, 3), 'markersize', 10), hold on, 
    yyaxis right
    plot(2, sg.metrics.r_asymp(k), 'ko', 'markerfacecolor', grayOrder(k)*ones(1, 3), 'markersize', 10),
end
yyaxis left, ylabel('time (ms)'), ylim([0, 180]), yyaxis right, ylim([0, 0.3])
xlim([0.5, 2.5]),  axis square, 
set(gca, 'xtick', [1, 2], 'xticklabel', {'t2pk', 't-asymp'}, 'fontsize', 14), box off

%% PART 2: FIT DN MODEL TO MUA DATA

m_mua = [];
init    = [0.05, 0.1, 1.8, 0.03, 0.01]; % 'tau1',  'tau2', 'n', 'sigma', 'shift'
lb      = [0.01,0.03,0,0,0];
ub      = [1, 1, 10, 1, 1];

t_mua = linspace(1/srate, size(dt_mua{1}, 2)./srate, size(dt_mua{1}, 2));

stim_mua = zeros(1, length(t_mua)); % stimulus on
stim_mua(400 : 900) = 1;
ext_stim_mua = [stim_mua, zeros(1, 300)];
ext_t_mua    = linspace(1/srate, length(ext_stim_mua)./srate, length(ext_stim_mua));


mua = struct();

for k = 1 : 2
    m_mua(k, :) = normMax_range(mean(dt_mua{k}), 100 : 300);
    % fit model
    mua.param(k, :) = fminsearchbnd(@(x) dn_computeFineFit(x, m_mua(k, :), stim_mua, t_mua, 'uniphasic'), init, lb, ub);
    mua.prd(k, :)   = normMax_range(dn_DNmodel([mua.param(k, 1), 0, mua.param(k, 2 : end), 1], ext_stim_mua, ext_t_mua), [100 : 300]);
end

mua.metrics = dn_computeDerivedParams(mua.param, 'uniphasic');

%% PART 2: VISUALIZATION

figure (1), 
for k = 1 : 2
   subplot_tight(4, 4, k + 4, plot_gap), cla, 
   patch([0.4, 0.9, 0.9, 0.4], [0, 0, 1, 1], 'k', 'FaceAlpha', 0.2, 'edgecolor', 'w'), hold on
   stem(t_mua, m_mua(k, :), 'k-', 'markersize', 1), 
   plot(ext_t_mua, mua.prd(k, :), 'r-', 'linewidth', 2)
   box off, set(gca, 'fontsize', 14), axis square, xlim([0.2, 1.2]), title('MUA')
   set(gca, 'xtick', [0.4, 0.9], 'xticklabel', [0, 0.5],  'ytick', [0, 1]), ax = gca; ax.XAxisLocation = 'origin'; ylim([-0.2, 1])
end

subplot_tight(4, 4, 8, plot_gap), 
for k = 1 : 2
    yyaxis left
    plot(1, mua.metrics.t2pk(k), 'ko', 'markerfacecolor', grayOrder(k)*ones(1, 3), 'markersize', 10), hold on, 
    yyaxis right
    plot(2, mua.metrics.r_asymp(k), 'ko', 'markerfacecolor', grayOrder(k)*ones(1, 3), 'markersize', 10),
end
yyaxis left, ylabel('time (ms)'), ylim([0, 180]), yyaxis right, ylim([0, 0.3])
xlim([0.5, 2.5]),  axis square, 
set(gca, 'xtick', [1, 2], 'xticklabel', {'t2pk', 't-asymp'}, 'fontsize', 14), box off

%% PART 3: FIT DN MODEL TO LFP BROADBAND DATA

lfp = struct();

% extract broadband form the lfp data
bands     = {[100, 170], 20}; % {[lower bound,  upper bound], window sz}
lineNoise =  50;
range     = [200 : 850];
ext_range = 200 : length(ext_t_mua);

for k = 1 : 2
    % notch filter
    lfp.dt{k} = dn_ecogNotch(dt_lfp{k}', srate, lineNoise);
    % extract broadband
    lfp.bb{k} = extractBroadband(lfp.dt{k}, srate, 4, bands);
    lfp.mbb(k, :) = normMax_range(median(lfp.bb{k}, 2), 200 : 400);
    
    % fit model
    lfp.param(k, :) = fminsearchbnd(@(x) dn_computeFineFit(x, lfp.mbb(k, range), stim_mua(range), t_mua(range), 'uniphasic'), init, lb, ub);
    lfp.prd(k, :)   = normMax_range(dn_DNmodel([lfp.param(k, 1), 0, lfp.param(k, 2 : end), 1], ext_stim_mua(ext_range), ext_t_mua(ext_range)), [1 : 200]);
end

lfp.metrics = dn_computeDerivedParams(lfp.param, 'uniphasic');

%% PART 3: VISUALIZE

figure (1), 
for k = 1 : 2
   subplot_tight(4, 4, k + 8, plot_gap), cla
   patch([0.4, 0.9, 0.9, 0.4], [0, 0, 1, 1], 'k', 'FaceAlpha', 0.2, 'edgecolor', 'w'), hold on
   plot(t_mua(range), lfp.mbb(k, range), 'b-', 'linewidth', 2)
   plot(ext_t_mua(ext_range), lfp.prd(k, :), 'r-', 'linewidth', 2)
   axis square, ax = gca; ax.XAxisLocation = 'origin'; ylim([-0.2, 1]), box off, set(gca, 'fontsize', 14)
   set(gca, 'xtick', [0.4, 0.9], 'xticklabel', [0, 0.5], 'ytick', [0, 1]), xlim([0.2, 1.2]), title('LFP broadband')
end

subplot_tight(4, 4, 12, plot_gap), 
for k = 1 : 2
    yyaxis left
    plot(1, lfp.metrics.t2pk(k), 'ko', 'markerfacecolor', grayOrder(k)*ones(1, 3), 'markersize', 10), hold on, 
    yyaxis right
    plot(2, lfp.metrics.r_asymp(k), 'ko', 'markerfacecolor', grayOrder(k)*ones(1, 3), 'markersize', 10),
end
yyaxis left, ylabel('time (ms)'), ylim([0, 180]), yyaxis right, ylim([0, 0.3])
xlim([0.5, 2.5]),  axis square, 
set(gca, 'xtick', [1, 2], 'xticklabel', {'t2pk', 't-asymp'}, 'fontsize', 14), box off

%% PART 4: FIT DN MODEL TO ECOG BROADBAND DATA
ecog = [];

t_ecog = [1 : size(dt_ecog, 3)]./1000;
stim_ecog = zeros(1, length(t_ecog)); stim_ecog(201 : 700) = 1;

for k = 1 : 2
   ecog.mbb(k, :) = normMax_range(squeeze(median(dt_ecog(k, :, :), 2)), [1 : 200]);
    % fit model
    ecog.param(k, :) = fminsearchbnd(@(x) dn_computeFineFit(x, ecog.mbb(k, :), stim_ecog, t_ecog, 'uniphasic'), init, lb, ub);
    ecog.prd(k, :)   = normMax_range(dn_DNmodel([ecog.param(k, 1), 0, ecog.param(k, 2 : end), 1], stim_ecog, t_ecog), [1 : 200]);
end

ecog.metrics = dn_computeDerivedParams(ecog.param, 'uniphasic');
%% PART 4: VISUALIZE

figure (1), 
for k = 1 : 2
    subplot_tight(4, 4, k + 12, plot_gap), cla
    patch([0.2, 0.7, 0.7, 0.2], [0, 0, 1, 1], 'k', 'FaceAlpha', 0.2, 'edgecolor', 'w'), hold on
    plot(t_ecog, ecog.mbb(k, :), 'b-', 'linewidth', 2)
    plot(t_ecog, ecog.prd(k, :), 'r-', 'linewidth', 2)
    axis square, xlim([0, 1]), ax = gca; ax.XAxisLocation = 'origin'; ylim([-0.2, 1]), box off, set(gca, 'fontsize', 14)
    set(gca, 'xtick', [0.2, 0.7], 'xticklabel', [0, 0.5], 'ytick', [0, 1]), title('ECoG broadband')
end

subplot_tight(4, 4, 16, plot_gap), 
for k = 1 : 2
    yyaxis left
    plot(1, ecog.metrics.t2pk(k), 'ko', 'markerfacecolor', grayOrder(k)*ones(1, 3), 'markersize', 10), hold on, 
    yyaxis right
    plot(2, ecog.metrics.r_asymp(k), 'ko', 'markerfacecolor', grayOrder(k)*ones(1, 3), 'markersize', 10),
end
yyaxis left, ylabel('time (ms)'), ylim([0, 180]), yyaxis right, ylim([0, 0.3])
xlim([0.5, 2.5]),  axis square, 
set(gca, 'xtick', [1, 2], 'xticklabel', {'t2pk', 't-asymp'}, 'fontsize', 14), box off

%% PART 5: BRIDGING SINGLE UNIT AND MUA

sg.nCells  = size(dt_sg, 2);
sg.cellNum = [1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 3];
nBoots     = 100;
sg.dt_boot = [];

% bootstrap over cells - extend data
dt_extend_sg = dt_sg(:, 1);

for k = 2 : sg.nCells
    tmp = repmat(dt_sg(:, k), [1, sg.cellNum(k)]); dt_extend_sg = [dt_extend_sg, tmp];
end

% bootstrap over cells - bootstrap and fit model
for k = 1 : nBoots
    idx = randi(sum(sg.cellNum), [1, sum(sg. cellNum)]); 
    sg.dt_boot(k, :) = median(dt_extend_sg(:, idx), 2);
    
    % fit model
    sg.groupParam(k, :) = fminsearchbnd(@(x) dn_computeFineFit(x, sg.dt_boot(k, :), stim_sg, t_sg, 'uniphasic'), init, lb, ub);
    sg.groupPrd(k, :)   =  normMax_range(dn_DNmodel([sg.groupParam(k, 1), 0, sg.groupParam(k, 2 : end), 1], stim_ecog, t_ecog), [1 : 200]);
end

%%

% make single unit prediction to ECoG stimuli
sg.groupMParam = median(sg.groupParam);
sg.groupMPrd   =  normMax_range(dn_DNmodel([sg.groupMParam(1), 0, sg.groupMParam(2 : end), 1], stim_ecog, t_ecog), [1 : 200]);

% make MUA prediction to ECoG stimuli
mua.mparam = median(mua.param);
mua.muaPrd = normMax_range(dn_DNmodel([mua.mparam(1), 0, mua.mparam(2 : end), 1], stim_ecog, t_ecog), [1 : 200]);

% make LFP prediction to ECoG stimuli
lfp.mparam = mean(lfp.param);
lfp.ecogPrd = normMax_range(dn_DNmodel([lfp.mparam(1), 0, lfp.mparam(2 : end), 1], stim_ecog, t_ecog), [1 : 200]);

%% PART 5: VISUALIZATION 

figure (2), clf
subplot_tight(3, 2, 1, plot_gap), 
plot(t_sg, normMax_range(dt_sg(:, 11), 1 : 20), 'k'), hold on
plot(t_sg, normMax_range(median(sg.dt_boot), [1 : 20])), axis square

subplot_tight(3, 2, 2, plot_gap)
plot(t_ecog, sg.groupMPrd, 'r-', 'linewidth', 2), hold on
plot(t_ecog, mua.muaPrd, 'k-', 'linewidth', 2)
axis tight, axis square, legend('single unit prediction', 'MUA')

subplot_tight(3, 2, 4, plot_gap), 
plot(t_ecog, lfp.ecogPrd, 'r-', 'linewidth', 2), hold on
plot(t_ecog, normMax_range(mean(ecog.prd), 1 : 200), 'k-', 'linewidth', 2)
legend('LFP prediction', 'ECoG')

axis tight, axis square

%% PART 6: BRIDGING MUA AND LFP broadband

%% PART 7: BRIDGING LFP broadband and ECoG broadband