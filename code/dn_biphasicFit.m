function target = trf_dcts_biphasicFit(param, data, stim, t)

normMax = @(x) x./max(x);

% compute model response
extprm = [param(1), param(2), param(3), 2, param(4), param(5), 1];

pred   = normMax(trf_dCTSmodel(extprm, stim, t));

target = sum((pred - data).^2);
% 
% figure (1), clf
% plot(pred), hold on
% plot(data),
% pause (0.5), drawnow


end