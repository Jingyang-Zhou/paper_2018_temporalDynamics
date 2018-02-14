% dn_chooseElectrodes

%% PREDEFINED VARIABLES

prm = [];

prm.ecog.roiNm = {'V1', 'V2', 'V3', 'lateral', 'ventral', 'dorsal'};
prm.ecog.selThresh = 2;
prm.ecog.tBase     = 1 : 200;

%% USEFUL FUNCTIONS

normBase = @(x, range) x - mean(x(range));

%% LOAD DATA

dataLoc = fullfile(dn_ECoG_RootPath, 'data');
fName   = 'dn_preprocessedData.mat';

a    = load(fullfile(dataLoc, fName));
data = a.dt.ecog;

%% DERIVED VARIABLES

label  = data.labels;
elec   = data.chans;
bb     = data.bb;

t_lth  = size(bb{1}, 2);
nElec  = size(bb{1}, 3);

prm.ecog.srate = 1000;
prm.ecog.t     = 1/prm.ecog.srate : 1/prm.ecog.srate : t_lth/prm.ecog.srate;

roiNms = {'V1', 'V2', 'V3', 'V3ab', 'hV4', 'VO', 'LO', 'TO', 'IPS'};
roiIdx = {1, 2, 3, 7, [5, 6], [4, 8]};

%% CHANGE TIME COURSE BROADBAND SIGNAL TO PERCENTAGE CHANGE

mbb{1} = squeeze(mean(bb{1}));
mbb{2} = squeeze(mean(bb{2}));

%% MODIFY ELECTRODE NAMES 

% Hand-modify the electrodes the names of which do not match the cortical
% locations. (VO and hV4 labels)
chans_mod = [103, 104, 99, 100, 94, 95, 57, 63,  62,  59, 105, 117];

for k = 1: length(chans_mod)
   modidx(k) = find(elec == chans_mod(k)); 
end

labs_mod = {'VO', 'VO', 'hV4', 'hV4', 'hV4', 'hV4', 'V3', 'LO', 'LO', 'LO',  '', ''};
label(modidx) = labs_mod;

%% CHANGE TIME COURSES TO PERCENTAGE SIGNAL

for iElec = 1 : nElec
    for ibands = 1 : 2
        % CONVERT THE SIGNAL TO PERCENTAGE CHANGE FROM THE BASELINE
        baseline = mean(mbb{ibands}(prm.ecog.tBase, iElec));
        mbb{ibands}(:, iElec) = mbb{ibands}(:, iElec)./baseline;
        % SUBTRACT THE BASELINE
        mbb{ibands}(:, iElec) = mbb{ibands}(:, iElec) - baseline;
    end
end

%% SELECT ELECTRODES BASED ON A THRESHOLD CRITERION






