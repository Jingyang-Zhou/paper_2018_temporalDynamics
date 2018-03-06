function scale = dn_fitScale(data, prd)
% DESCRIPTION -------------------------------------------------------------
% function scale = dn_fitScale(data, prm, stim, t, irfType)
%
% INPUTS ------------------------------------------------------------------
% param:
% data: nBoots x time course
% stim:
% t: 
% irfType: 'uniphasic' and 'biphasic'
%
% OUTPUTS -----------------------------------------------------------------
% scale:

%% USEFUL FUNCTIONS

normMax = @(x) x./max(x);

%% COMPUTE SCALE

for k = 1 : size(data, 1)  
    % NORMALIZE THE PREDICTIONS
    prd(k, :)  = normMax(prd(k, :));
    
    % COMPUTE THE SCALE
    scale(k) = regress(data(k, :)', prd(k, :)');
end

end