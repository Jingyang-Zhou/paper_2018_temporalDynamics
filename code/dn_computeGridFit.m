function r2 = dn_computeGridFit(prms, data, stim, t, irfType)
% function r2 = dn_computeGridFit(prms, data, stim, t, irfType):
% This functions handles fitting the DN model.
%
% INPUTS ------------------------------------------------------------------
% prms    : DN model parameters. If the irfType is "uniphasic", then "prms" 
%           contains 4 model parameters in the order of tau1, tau2, n and sigma. 
%           If the irfType is "biphasic", then four different model
%           parameters are expected: tau1, weight, tau2, sigma.
% data    : nTrials(nBoots) x time course
% stim    : stimulus, a time course with 1 represents contrast increment
%           and 0 otherwise.
% t       : time, in unit of seconf
% irfType : can either be "uniphasic" or "biphasic".
%
% OUTPUTS -----------------------------------------------------------------
% r2      : the squared correlation between the model prediction and data

%% PRE-DEFINE PARAMETERS/FUNCTIONS

normMax = @(x) x./max(x);

%% compute model response

% For the dn_DNmodel input parameters: {'tau1', 'weight', 'tau2', 'n', 'sigma', 'shift', 'scale'};

% SET UP INPUT MODEL PARAMETERS -------------------------------------------
switch irfType
    case 'uniphasic'
        modelPrms = [prms(1), 0, prms(2 : 4)', 0, 1];
    case 'biphasic'
        modelPrms = [prms(1 : 3)', 2, prms(4), 0.05, 1];
end

% COMPUTE MODEL PREDICTIONS -----------------------------------------------
rsp = normMax(dn_DNmodel(modelPrms, stim, t));

% COMPUTE R2^2 BETWEEN MODEL PREDICTIONS AND DATA -------------------------
r2 = corr(rsp', data').^2;

%% visualize

% figure (1), clf
% plot(normMax(rsp), 'r-'), hold on
% plot(normMax(data), 'b-'), drawnow, pause(0.1)

end