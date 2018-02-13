% compare single unit, mua and broadband signals

%%  generate single unit predictions

stim = [zeros(1, 200), ones(1, 500), zeros(1, 500)];
t    = [1 : length(stim)]./1000;

norm_max = @(x) x./max(x);

% tau1, weight, tau2, n, sigma, shift, scale
sparam = [0.07, 0, 0.07, 2.8, 0.014, 0.03, 1];

original = norm_max(dn_DNmodel(sparam, stim, t));

%% perturb parameters

generated =[];

nsamples = 20;

for k = 1 : nsamples
   r_param(k, :) = sparam +[max(0, randn*0.2), 0, max(0, randn*0.2), max(0, randn*0.2), max(0, randn*0.02), 0, 0]; 
   generated(k, :) = norm_max(dn_DNmodel(r_param(k, :), stim, t));
end

% plot
figure (1), %clf
subplot(1, 2, 1)
plot(t, stim, 'k:', t, original, 'b-'), hold on
plot(t, norm_max(mean(generated)), 'r-')
subplot(1, 2, 2)
plot(t, generated)

%% perturb parameters 2

clear generated

nsizes = [1, 1000];

for k0 = 1 : length(nsizes)
    r_param = {};
    for k = 1 : nsizes(k0)
        r_param{1}(k, :) = abs(sparam + [randn*sparam(1)*1.5, 0, 0, 0, 0, 0, 0]); % tau1
        r_param{2}(k, :) = abs(sparam +[0, 0, 0, 0, 0, randn*sparam(6)*1.5, 0]); % tau2
        r_param{3}(k, :) = abs(sparam +[randn*sparam(1), 0, 0, randn*sparam(4)*0.5, randn*sparam(5), randn*sparam(6), 0]); % n
        for k1= 1 : 3 % which parameter is distributed
            generated{k0}(k1, k, :) =  norm_max(dn_DNmodel(r_param{k1}(k, :), stim, t));
        end
    end
    m_generated{k0} = squeeze(mean(generated{k0}, 2));
    m_generated{k0} = norm_max(m_generated{k0}');
end
%
figure (2), clf
subplot(3, 1, 1), plot(t, stim, 'k:', t, original, 'k-'),
for k = 1 : 2
    subplot(3, 1, k + 1)
    plot(t, stim, 'k:', t, original, 'k-'), hold on

    plot(t, m_generated{k}), box off
  %  title(num2str(nsizes(k)))
  legend('stimulus', 'original', 'summation window', 'onset time', 'extent of normalization')
end

%% analyze the effect of perturbing parameters
lb = [0.01, 0.01, 0.5, 0.001, 0];
ub = [1, 1, 5, 1, 1];
sparam1 = [0.07, 0.07, 2.8, 0.014, 0.03];

nsamples2 = 500;

for k = 1 : nsamples2
    param_new{1}(k, :) = sparam + [max(0.01, randn *0.06), 0, 0, 0, 0, 0, 0]; 
    param_new{2}(k, :) = sparam + [0, 0,max(0.001, randn *0.069), 0, 0, 0, 0]; 
    param_new{3}(k, :) = sparam + [0, 0, 0, max(0.5, randn *1), 0, 0, 0];
    param_new{4}(k, :) = sparam + [0, 0, 0, 0, max(0.0001, randn *0.02), 0, 0];
end

% fit model
parfor k1 = 1 : 4
    for k = 1 : nsamples2
        % generate new data
        data = norm_max(trf_dCTSmodel(param_new{k1}(k, :), stim, t));
        param_fit{k1}(k, :) = fminsearchbnd(@(x) trf_dcts_finefit(x, data, stim, t), sparam1, lb, ub);
    end
    k1
end


%% plot parameters

figure (2), clf
subplot(2, 2, 1)
plot(param_new{1}(:, 1), param_fit{1}(:, 1), 'k.', 'markersize', 5), hold on
subplot(2, 2, 2)
plot(param_new{2}(:, 3), param_fit{2}(:, 2), 'k.', 'markersize', 5)
subplot(2, 2, 3)
plot(param_new{3}(:, 4), param_fit{3}(:, 3), 'k.', 'markersize', 5)
subplot(2, 2, 4)
plot(param_new{4}(:, 5), param_fit{4}(:, 4), 'k.', 'markersize', 5)