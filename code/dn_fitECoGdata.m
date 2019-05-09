% dn_fit ECoG data

%% pre-defined parameters

t_preStim = 1 : 200;

%% load data and parameters

dataLoc = fullfile(dn_ECoG_RootPath, 'data');
fName   = 'dn_data.mat';
a       = load(fullfile(dataLoc, fName));


nrois   = length(a.param.roiNm);
nBoots  = 100;

data  = reshape(permute(a.bsData, [2, 1, 3]), nrois * nBoots, []);
param   = a.param;
raw     = a.raw;
bsData  = a.bsData;

% make stimulus
stim = zeros(1, size(data, 2));
stim(201 : 700) = 1;
t = linspace(0.001, 1.204, size(data, 2));

if ~isfield(1, 'dn')
    dn = [];
else
    dn = a.dn;
end

% for debugging:
% data = data(1 : 10, :);

%% grid fit

irfType = 'uniphasic';

[dn.modelSeed, dn.seedR2] = dn_gridFit(data, param, stim, t, irfType);

%% fine fit
% temporary 
dn = a.dn;

seed = [dn.modelSeed, zeros(size(dn.modelSeed, 1), 1)];
[dn.param, dn.prd, dn.r2, dn.exitflg] = dn_fineFit(data, stim, t, param, seed, irfType);

%% fit the scale of the ECoG broadband response

scale    = dn_fitScale(data, dn.param, stim, t, irfType);
dn.param = [dn.param, scale'];

%% compute derived parameters

dn.derivedPrm = dn_computeDerivedParams(dn.param, irfType);

%% save dn

%% save file

dn_biphasic = a.dn_biphasic;

save(fullfile(dataLoc, fName), 'bsData', 'param', 'raw', 'dn', 'dn_biphasic')