function [diff, pred_ds] = trf_fitModel_contrast(params, data, time, contrast_levels, stim)
%function diff = trf_fitModel_contrast(params, data, time, contrast_levels, stim)

% params has n entries, tau1, tau2, n, sigma, shift

% make predictions using the parameters
extended_prm = [params(1), 0, params(2), params(3), params(4), params(5), 1];
% time = time - time(1) + 1;
 t = [1 : length(stim)]./1000;

for k = 1 : length(contrast_levels)
   pred(k, :) = dn_DNmodel(extended_prm, stim.*contrast_levels(k), t);
end
% down_sample predictions
pred_ds = pred(:, time);

normMax = @(x) x./max(x(:));
pred_ds = normMax(pred_ds);

% compare to data
diff = sum((pred_ds(:) - data(:)).^2);


%% visualize
% figure (100), clf
% 
% plot(time, data, 'r-'), hold on
% plot(time, pred_ds, 'k-'), drawnow

end