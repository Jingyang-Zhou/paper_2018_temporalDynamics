function [normrsp, linrsp, numrsp, demrsp] = dn_DNmodel(param, stim, t)
%
% function normrsp = dn_DNmodel(param, stim, t)
% INPUTS :
% params : 6 fields.
%          tau1 -- irf peak time, in unit of second
%          weight -- the weight in the biphasic irf function, set weight to
%          0 if want to use uniphasic irf function.
%          tau2 -- time window of adaptation, in unit of second
%          n -- exponent
%          sigma -- semi-saturation constant
%          shift -- time between stimulus onset and when the signal reaches
%          the cortex, in unit of second
%          scale -- response gain.
% OUTPUTS:
% normrsp: normalization response.

% 04/05 I added parameter field "n", 

%% pre-defined variables

x   = [];
%x.n = 2;

normSum = @(x) x./sum(x);

%% initiate model fitting

fields = {'tau1', 'weight', 'tau2', 'n', 'sigma', 'shift', 'scale'};
x      = toSetField(x, fields, param);

%% compute response

% compute irf
irf1 = gammaPDF(t, x.tau1, 2);
irf2 = gammaPDF(t, x.tau1*1.5, 2);
irf  = irf1 - x.weight.* irf2;

% compute pool filter
irf_norm = normSum(exp(-t/x.tau2));

ctLth = length(irf);


for istim = 1 : size(stim, 1)
    % modify stimulus based on shift
    dt        = t(2) - t(1);
    sft       = round(x.shift * (1/dt));
    stimtmp   = padarray(stim(istim, :), [0, sft], 0, 'pre');
    stim(istim, :) = stimtmp(1 : size(stim, 2));
    
    % compute linear response
    linrsp(istim, :)  = convCut(stim(istim, :), irf, ctLth);
    % compute normalization pool response
    poolrsp(istim, :) = convCut(linrsp(istim, :), irf_norm, ctLth);
    % compute numerator:
    numrsp(istim, :)  = linrsp(istim, :).^x.n;
    % compute denominator:
    demrsp(istim, :)  = x.sigma.^x.n + poolrsp(istim, :).^x.n;
    
    % compute normalization response
    normrsp(istim, :) = x.scale.*(numrsp(istim, :)./demrsp(istim, :));
end

end