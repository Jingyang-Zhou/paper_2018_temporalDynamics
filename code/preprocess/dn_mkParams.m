%% dn_mkParam

%% make parameteres

param = [];

param.exp.indiRoiLabel = {'V1', 'V2', 'V3', 'LO', 'hV4', 'VO', 'V3A', 'IPS'};
param.exp.cmbRoiLabel  = {'V1', 'V2', 'V3', 'LA', 'VA', 'DA'}; % lateral, ventral and dorsal areas
param.exp.stimType     = {'w', 'p', 'b'}; % for white, pink and brown noise

%% relative to model fitting

param.model.nFit     = 1000;
param.model.tau1rnd  = 0.01 + (1 - 0.01).*rand(1, nFit);   % between 0.01 and 1
param.model.tau2rnd  = 0.01 + (1 - 0.01).*rand(1, nFit);   % between 0.01 and 1
param.model.nrnd     = 0.5 + (5 - 0.5).*rand(1, nFit);     % between 0.1 and 5
param.model.sigmarnd = 0.01 + (0.5 - 0.01).*rand(1, nFit); % between 0.01 and 0.5
param.model.lower    = [0.07, 0.01, 0, 0.01, 0.001]; % lower bound for the model parameters
param.model.upper    = [1, 2, 6, 1, 0.2]; % upperbound for the model parameters
prm.model.fitOptions = optimset('Display', 'notify', 'Algorithm', 'trust-region-reflective', 'MaxIter', 5000);
