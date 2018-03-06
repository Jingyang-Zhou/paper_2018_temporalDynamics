function dif = dn_fitExponential(tau, data, range, t)


%% MAKE MODEL PREDICTION

pred = exp(-t./tau);

%% COMPUTE THE DIFFERENCE BETWEEN DATA AND MODEL PREDICTION

rng_data = data(range);
rng_pred = normMax(pred(range));

dif = sum((rng_data - rng_pred).^2);

%% VISUALIZE THE MODEL FIT

% figure (100), clf
% plot(rng_data), hold on
% plot(rng_pred), drawnow

end