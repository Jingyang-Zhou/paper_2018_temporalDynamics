% dn_makeFigureS6

%% useful functions

normMax_range = @(x, range) (x - mean(x(range)))./max(x - mean(x(range)));
normMax       = @(x) x./max(x);

%% load single unit data

dataLoc = fullfile(dn_ECoG_RootPath, 'data');
fName   = 'Figure4Data.xlsx';

a = xlsread(fullfile(dataLoc, fName));

t      = a(:, 1)./1000;
stim   = ones(1, length(t));

data   = a(:, 2 : end);

nCells = size(data, 2);

cellNum = [1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 3];

whichcells = 1 : 12;

%% load ecog data

dataLoc = fullfile(dn_ECoG_RootPath, 'data');
fName   = 'dn_data.mat';

a       = load(fullfile(dataLoc, fName));
param   = a.param;
dn      = a.dn;

t_ecog    = [1 : 1200]./1000;
stim_ecog = zeros(1, length(t_ecog)); stim_ecog(201 : 400) = 1;

%% compare 3 most common individual cell response types and ECoG summary metrics

figure (1), clf
for k = 1 : nCells
    subplot(4, 3, k)
    plot(t, data(:, k), 'k-'), box off, title(cellNum(k))
    set(gca, 'ytick', [0, 1]), axis square
end

%% compare predicted MUA and ECoG summary metrics
repCell = 1 : 12;
repData = data(:, repCell);
sua     = [];

% fit model to the 3 representative cells:
%init = [0.1, 0.1, 2, 0.1, 0];
lb   = [0.03,0.03,0,0.01,0];
ub   = [1, 1, 6, 0.5, 0.1];

% STEP (1): GRID SEARCH
for k = 1 : length(repCell)
    [sua.REPseed(k, :), sua.REPseedR2(k)] = dn_gridFit(repData(:, k)', param, stim, t, 'uniphasic');
end

% STEP (2): FINE FIT
for k = 1 : length(repCell)
    [sua.REPparam(k, :), ~, exitflg(k)] = fminsearchbnd(@(x) dn_computeFineFit(x, repData(:, k)', stim, t, 'uniphasic'), [sua.REPseed(k, :), 0], lb, ub);
    sua.REPpred(k, :)  = normMax(dn_DNmodel([sua.REPparam(k, 1), 0, sua.REPparam(k, 2 : end), 1], stim, t));
    sua.REPECoGPred(k, :) = normMax(dn_DNmodel([sua.REPparam(k, 1), 0, sua.REPparam(k, 2 : end), 1], stim_ecog, t_ecog));
end

% STEP (3): COMPUTE SUMMARY METRICS
sua.REPsummary = dn_computeDerivedParams(sua.REPparam, 'uniphasic');

%% Visualize individual cell and ECoG response and summary metrics

% plot summary metrics (compare single cell and ECoG response)
ecog_mt2pk  = median(dn.derivedPrm.t2pk(1 : 100));
ecog_masymp = median(dn.derivedPrm.r_asymp(1 : 100));

ecog_st2pk  = prctile(dn.derivedPrm.t2pk(1 : 100), [25, 75]);
ecog_sasymp = prctile(dn.derivedPrm.r_asymp(1 : 100), [25, 75]);

% single unit parameters
sua_mt2pk  = median(sua.REPsummary.t2pk);
sua_masymp = median(sua.REPsummary.r_asymp);

sua_st2pk  = prctile(sua.REPsummary.t2pk, [25, 75]);
sua_sasymp = prctile(sua.REPsummary.r_asymp, [25, 75]);


figure (2), clf
% ecog parameters
yyaxis left, plot(1, ecog_mt2pk, 'ko', 'markersize', 12, 'markerfacecolor', 'k'), hold on, plot([1, 1], ecog_st2pk, 'k-')
yyaxis right, plot(2, ecog_masymp, 'ko', 'markersize', 12, 'markerfacecolor', 'k'), hold on, plot([2, 2], ecog_sasymp, 'k-')
% single unit parameters
yyaxis left, plot(1.2, sua.REPsummary.t2pk, 'ro')
yyaxis right, plot(2.2, sua.REPsummary.r_asymp, 'ro')
%yyaxis left, plot(1.2, sua_mt2pk , 'ro', 'markersize', 12, 'markerfacecolor', 'r'), hold on, plot([1.2, 1.2], sua_st2pk , 'r-')
%yyaxis right, plot(2.2, sua_masymp, 'ro', 'markersize', 12, 'markerfacecolor', 'r'), hold on, plot([2.2, 2.2], sua_sasymp, 'r-')

xlim([0.5, 2.5]), %yyaxis left, ylim([0, 240]), yyaxis right, ylim([0, 1])

