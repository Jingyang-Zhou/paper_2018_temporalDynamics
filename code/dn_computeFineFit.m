function [r2, rsp] = dn_computeFineFit(prms, data, stim, t, irfType)
% DESCRIPTION -------------------------------------------------------------
%
% INPUTS ------------------------------------------------------------------
% prms    : DN model parameters. If the irfType is "uniphasic", then "prms" 
%           contains 5 model parameters in the order of tau1, tau2, n sigma, and shift. 
%           If the irfType is "biphasic", then 5 different model
%           parameters are expected: tau1, weight, tau2, sigma, and shift.
% data    : 1 x time course
% stim    : stimulus, a time course with 1 represents contrast increment
%           and 0 otherwise.
% t       : time, in unit of seconf
% irfType : can either be "uniphasic" or "biphasic"
%
% OUTPUTS -----------------------------------------------------------------
% r2      : sum of the squared difference between data and model prediction

%% USEFUL FUNCTION

normMax = @(x) x./max(x);

%% PRE-DEFINED VARIABLE

data = normMax(data);

%% COMPUTE MODEL PREDICTIONS

prms = dn_fillDNParams(prms, irfType);
rsp  = normMax(dn_DNmodel(prms, stim, t));

%% COMPUTE R2

r2 = sum((data - rsp).^2);

%% VISUALIZE

% figure (1), clf
% plot(rsp, 'r-'), hold on
% plot(data, 'b-'), drawnow

end

