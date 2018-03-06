function [seed, seedR2] = dn_gridFit(data, param, stim, t, irfType)
% DESCRIPTION -------------------------------------------------------------
% function [seed, seedR2] = dn_gridFit(data, param, stim, t, irfType)
% This file compute the DN model prediction over a parameter grid, and then
% it extracts the set of the parameters on the grid that make the best
% prediction and suggests the parameter set as seed for the fine model fit.
%
% INPUTS ------------------------------------------------------------------
% data   : n bootstraps x time course
% irfType: can either be "uniphasic" or "biphasic
% param  : a struct that contains the parameter steps to make the parameter
%          grid.
%
% OUTPUT ------------------------------------------------------------------
% seed   : the set of parameters that we use as a seed for the fine fit
% seedR2 : the R2 of the model seed predictions 
% 
% DEPENDENCIES ------------------------------------------------------------
% dn_computeGridFit.m
% dn_DNmodel.m

%% MAKE GRID

switch irfType
    case 'uniphasic'
        tau1 = param.tau1steps; tau2 = param.tau2steps; n = param.nsteps; sigma = param.sigmasteps;
        [s{1}, s{2}, s{3}, s{4}] = ndgrid(tau1, tau2, n, sigma);
    case 'biphasic'
        s{1} = param.weightsteps;
end

%% DERIVED PARAMETERS

nBoots = size(data, 1);
nSteps = length(s{1}); % Here assumes that each parameter has equal number of steps

%% RESHAPE THE GRID TO A LONG VECTOR

nPrm = length(s); gridLth = nSteps^nPrm; grid = [];

% If each parameter has 10 steps, here is to make the grid parameter of
% size [1, 10^4], when there are 4 parameters.
switch irfType
    case 'uniphasic'
        for iPrm = 1 : nPrm
            grid(iPrm, :) = reshape(s{iPrm}, [1, gridLth]);
        end
    case 'biphasic'
        grid = s{1};
end

%% COMPUATE GRID PREDICTIONS
r2 = [];
switch irfType
    case 'uniphasic'
        for iBoot = 1 : nBoots
            for iSet = 1 : gridLth
                r2(iBoot, iSet) = dn_computeGridFit(grid(:, iSet), data(iBoot, :), stim, t, irfType);
            end
            iBoot
        end
    case 'biphasic'
        for iBoot = 1 : nBoots
            for iSet = 1 : gridLth
                thisprm = [0.02, s{1}(iSet), 0.1, 0.05]';
                r2(iBoot, iSet) = dn_computeGridFit(thisprm, data(iBoot, :), stim, t, irfType);
            end
        end
    otherwise, error('Unidentifiable irftype.')
end

%% FIND THE SET OF PARAMETERS THAT PRODUCE THE MAXIMUM R2

for iBoot = 1 : nBoots
    this_r2 = r2(iBoot, :);
    % FIND THE MAXIMUM R2 -------------------------------------------------
    seedR2(iBoot)     = max(this_r2);
    seed_r2Idx(iBoot) = find(this_r2 == seedR2(iBoot), 1);
    seed(iBoot, :)    = grid(:, seed_r2Idx(iBoot));
end

end