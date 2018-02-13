function [modelSeed, modelR2] = dn_gridFit(data, param, stim, t, irfType)
% DESCRIPTION
%
% INPUTS ---------------------------------------
%
% irfType: can either be "uniphasic" or "biphasic
% param  : experimental parameters
%
% OUTPUT ---------------------------------------
% 
% DEPENDENCIES ---------------------------------

%% make grids

% TEMPORARY CODE HERE, TO BE REMOVED

switch irfType
    case 'uniphasic'
        [s{1}, s{2}, s{3}, s{4}] = ndgrid(param.tau1steps, param.tau2steps, param.nsteps, param.sigmasteps);
    case 'biphasic'
        tau1steps = linspace(0.03, 0.15, 10);
        tau2steps = linspace(0.03, 0.15, 10);
        weightsteps = linspace(0, 1, 10);
        sigmasteps = linspace(0.001, 0.2, 10);
        % [s{1}, s{2}, s{3}, s{4}] = ndgrid(param.tau1steps, param.weightsteps, param.tau2steps, param.sigmasteps);
        [s{1}, s{2}, s{3}, s{4}] = ndgrid(tau1steps, weightsteps, tau2steps, sigmasteps);
end

% make the grid
grid = []; gridLth = length(s{1});

for k = 1 : length(s)
   grid(k, :) = reshape(s{k}, [1, gridLth^length(s)]); 
end

%% do grid fit

maxIdx    = [];
modelSeed = [];
modelR2   = [];

r2 = [];

for k = 1 : size(data, 1)
    for k1 = 1 : size(grid, 2)
        r2(k, k1) = dn_computeGridFit(grid(:, k1), data(k, :), stim, t, irfType);
    end
end

%% Find the set of parameters that produce the maximum r2 

for k = 1 : size(data, 1)
    maxIdx(k)       = find(r2(k, :) == max(r2(k, :)), 1);
    modelSeed(k, :) = grid(:, maxIdx(k));
    modelR2(k, :)   = r2(k, maxIdx(k));
end

end