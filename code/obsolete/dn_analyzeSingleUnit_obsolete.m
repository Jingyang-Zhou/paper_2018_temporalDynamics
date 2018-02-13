% fit single unit data

%% useful function

normMax_range = @(x, range) (x - mean(x(range)))./max(x - mean(x(range)));
normMax       = @(x) x./max(x);

%% load data 

dataLoc = fullfile(temporalRootPath, 'data', 'ECoG');
fName   = 'figure4Data.xlsx';

a = xlsread(fullfile(dataLoc, fName));

t     = a(:, 1);
data  = a(:, 2 : end);
mdata = normMax_range(median(data, 2), 1 : 20);

nCells = size(data, 2);

%% visualize the data

figure (1), clf

for k = 1 : nCells
   subplot(4, 3, k)
   plot(t, data(:, k), 'linewidth', 2), hold on
   plot(t, mdata, 'k:', 'linewidth', 2), box off
end

%% extend data and prep for bootstrap

% label number of cells in each panel
cellNum = [1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 3];

extend_data = data(:, 1);

for k = 2 : nCells
    tmp = repmat(data(:, k), [1, cellNum(k)]);
    extend_data = [extend_data, tmp];
end

%% bootstrap: 

nBoots = 100;

boot_mdata = [];

for k = 1 : nBoots
    idx = randi(sum(cellNum), [1, sum(cellNum)]);
    boot_mdata(k, :) = mean(extend_data(:, idx), 2);
end

figure (2), clf
plot(boot_mdata')

%% fit model to bootstrapped data

stim = ones(1, length(t));

init  = [0.04, 0.1, 1.8, 0.05, 0.01, 1]; % 'tau1',  'tau2', 'n', 'sigma', 'shift', 'scale'
lb    = [0,0,0,0,0,0];
ub    = [1, 1, 10, 1, 1, 10];

param = [];

for k = 1 : nBoots
    param(k, :) = fminsearchbnd(@(x) dn_computeFineFit(x, boot_mdata(k, :), stim, t./1000, 'uniphasic'), init, lb, ub);
end

%% check whether the averaged single unit spike rate is similar to mua

% load mua
dataLoc = '/Volumes/server/Projects/Temporal_integration/data/ECoG/';
dataNm  = 'CONTEXT.mat';

a = load(fullfile(dataLoc, dataNm));

% a.DETS: [electrode, session number, condition]
% data: 930 samples in length, 400-900 stimulus on

range = 2 : size(a.MUA, 2); % this is the range of data that we will be analyzing

% get electrode and session index
elecIdx = a.DETS(:, 1);
sessIdx = a.DETS(:, 2);

mua{1} = a.MUA(elecIdx == 6, range);
mua{2} = a.MUA(elecIdx == 7, range);

figure (1), clf
plot(normMax_range(median(mua{1}), [1 : 300])), hold on
plot(371 : 570, normMax(mean(boot_mdata)))
