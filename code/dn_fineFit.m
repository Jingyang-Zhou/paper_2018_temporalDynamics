function [modelprm, modelprd, modelr2, exitflg] = dn_fineFit(data, stim, t, param, seed, irfType)


%% prepare to fit model

switch irfType
    case 'uniphasic'
        %lb = param.lb;
        %ub = param.ub;
        lb = [0.07, 0.07, 0, 0.01, 0.001];
        ub = [1, 1, 6, 1, 0.2];
    case 'biphasic'
        %lb = param.bilb;
        %ub = param.biub;
        lb   = [0, 0, 0, 0.001, 0];
        ub   = [1, 1, 1, 1, 0.5];
end

%% fit model

modelprm = []; modelprd = []; modelr2 = [];

options = optimset('MaxFunEvals', 1000);

for k = 1 : size(data, 1)
    [modelprm(k, :), ~, exitflg(k)] = fminsearchbnd(@(x)dn_computeFineFit(x, data(k, :), stim, t, irfType), seed(k, :), lb, ub, options);
   
    % compute model prediction
    switch irfType
        case 'uniphasic'
            prms = [modelprm(k, 1), 0, modelprm(k, 2), modelprm(k, 3), modelprm(k, 4), modelprm(k, 5), 1];
        case 'biphasic'
            prms = [modelprm(k, 1), modelprm(k, 2), modelprm(k, 3), 2, modelprm(k, 4), modelprm(k, 5), 1];
    end
    modelprd(k, :) = dn_DNmodel(prms, stim, t);

    % compute model r2
    modelr2(k) = corr(modelprd(k, :)', data(k, :)').^2;
   
    % track progress
    if rem(k, 50) == 0, disp(sprintf('dn_FineFit: %d', k)), end
end

end