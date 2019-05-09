% Extract electrode eccentricity

%% LOAD PREPROCESSED DATA

dataLoc = fullfile(dn_ECoG_RootPath, 'data');
dataNm  = 'dn_preprocessedData.mat';

a = load(fullfile(dataLoc, dataNm));

%% DERIVED PARAMETERS

elec_roi = dt.ecog.elec_roi;
nrois    = length(elec_roi);

%% CHECK ELECTRODE ECCENTRICITY

b = load('/Volumes/server/Projects/ECoG/Figures/pRF_data/Power/oneBar/CAR/Subj17_exp_bb.mat');
x0 = b.params.params(:, 1);
y0 = b.params.params(:, 2);
s0 = b.params.params(:, 3);

subj = 17;

for iroi = 1 : nrois
    idx = elec_roi{iroi};
    [x{iroi}, y{iroi}, s{iroi}] = ecogFitPRFPix2Deg(subj, x0(idx), y0(idx), s0(idx));
    ecc{iroi} = sqrt(x{iroi}.^2 + y{iroi}.^2);
end

%% SAVE DATA

dt = a.dt;
dt.ecog.ecc = ecc;

save(fullfile(dataLoc, dataNm), 'dt')

 