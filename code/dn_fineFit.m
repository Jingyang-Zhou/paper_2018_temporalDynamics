function [prm, prd, r2, exitFlg] = dn_fineFit(data, stim, t, param, seed, irfType)
% DESCRIPTION -------------------------------------------------------------
% function [modelprm, modelprd, modelr2, exitflg] = dn_fineFit(data, stim, t, param, seed, irfType)
%
% INPUTS ------------------------------------------------------------------
% data : nRuns x time course
% stim : stimulus time course
% t    : time, in unit of second
% param: the structure that constains the upper and lower bound of the
%        uniphasic DN model fit
% seed : nRuns x seed values (uniphasic: tau1, tau2, n, sigma, shift)
% irfType : 'uniphasic' or 'biphasic'
%
% OUTPUTS -----------------------------------------------------------------
% prm
% prd
% r2
% exitFlg

%% PRE-DEFINED / EXTRACTED PARAMETERS

nBoots  = size(data, 1);
options = optimset('MaxFunEvals', 5000);

%% SET UP LOWER AND UPPWER BOUND FOR THE MODEL FIT

lb = param.lb; ub = param.ub;

%% FIT THE DN MODEL TO THE DATA

for iBoot = 1 : nBoots
    % FIT DN MODEL TO THE DATA --------------------------------------------
    thisdt = data(iBoot, :); thisSd = seed(iBoot, :);
    
    [prm(iBoot, :), ~, exitFlg(iBoot)] = ...
        fminsearchbnd(@(x)dn_computeFineFit(x, thisdt, stim, t, irfType), thisSd, lb, ub, options);
    % COMPUTE MODEL PREDICTION --------------------------------------------
    thisprm       = dn_fillDNParams(prm(iBoot, :), irfType);
    prd(iBoot, :) = dn_DNmodel(thisprm, stim, t);
    
    % COMPUTE MODEL R2 ----------------------------------------------------
    r2(iBoot) = corr(prd(iBoot, :)', data(iBoot, :)').^2;
    
   % TRACK FITTING PROGRESS -----------------------------------------------     
   if rem(iBoot, 50) == 0, disp(sprintf('dn_FineFit: %d', iBoot)), end
end

end