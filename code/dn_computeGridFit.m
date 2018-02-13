function r2 = dn_computeGridFit(prms, data, stim, t, irfType)

% irfType: can either be "uniphasic" or "biphasic"


% fields = {'tau1', 'weight', 'tau2', 'n', 'sigma', 'shift', 'scale'};

%% compute model response

switch irfType
    case 'uniphasic'
        modelPrms = [prms(1), 0, prms(2), prms(3), prms(4), 0, 1];
    case 'biphasic'
        modelPrms = [prms(1), prms(2), prms(3), 2, prms(4), 0, 1];
end

% compute model response
rsp = dn_DNmodel(modelPrms, stim, t);

% compute r2
r2 = corr(rsp', data').^2;


%% visualize

% figure (1), clf
% plot(rsp, 'r-'), hold on
% plot(data, 'b-'), drawnow

end