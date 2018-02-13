function [r2, rsp, rsp_normMax] = dn_computeFineFit(prms, data, stim, t, irfType)

% irfType: can either be "uniphasic" or "biphasic"


% fields = {'tau1', 'weight', 'tau2', 'n', 'sigma', 'shift', 'scale'};

%% compute model response

switch irfType
    case 'uniphasic'
        modelPrms = [prms(1), 0, prms(2), prms(3), prms(4), prms(5), 1];
    case 'biphasic'
        modelPrms = [prms(1), prms(2), prms(3), 2, prms(4), prms(5), 1];
end

% compute model response
normMax = @(x) x./max(x);

rsp         = dn_DNmodel(modelPrms, stim, t);
rsp_normMax = normMax(dn_DNmodel(modelPrms, stim, t));
data        = normMax(data);

% compute r2
%r2 = -corr(data', rsp_normMax');
r2 = sum((data - rsp_normMax).^2);

%% visualize
% % 
% figure (1), clf
% plot(rsp_normMax, 'r-'), hold on
% plot(data, 'b-'), drawnow

end

