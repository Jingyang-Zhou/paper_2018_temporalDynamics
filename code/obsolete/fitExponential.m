function dif = fitExponential(tau, multiUnit, lfp, t)


%% useful functions

normMax_range = @(x, range) (x - mean(x(range)))./max(x - mean(x(range)));

%% construct an exponential smoothing function

smth = exp(-t/tau);

%% make prediction

prd  = normMax_range(convCut(multiUnit, smth, length(t)), 1 : 50);
dif = sum((prd-lfp).^2)

%% visualize

% figure (1), clf
% plot(t, multiUnit, 'k-'), hold on
% plot(t, prd, 'b-')
% plot(t, lfp, 'r-'), drawnow  

end
