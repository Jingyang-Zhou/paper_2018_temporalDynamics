% fit to single and multi-unit data

% TO CONFIRM:

% In the Self data, the stimulus length is 200 ms!
% The meaning of the third and the fourth column of a.DETS

% Interpretation of the Self data:
% The first dimension of a.DETS represents which electrodes
% The second dimension of a.DETS represents a number from 1 to 121,
% representing the location of the stimulus
% The third dimension represents a random number from 4 to 14.
% The fourth dimension represents a random number from -8 to 2.

% pre-processing procedure:
% (1) separate the time seriese by electrodes
% (2) separate the time seriese by stimulus locations
% (3) threshold the time seriese


%% single unit: load data

fLoc  = fullfile(temporalRootPath, 'data', 'ECoG');
fName = 'RFmapping.mat';

a = load(fullfile(fLoc, fName));

%% pre-defined parameters

n_locations = 121;

criterion = 4;

norm_max = @(x) x./max(x);


%% pre-process data: normalize the time series

% use the 100 ms before the stimulus onset as a reference

mua = {};

t_idx = (a.TD >= -0.1) & (a.TD < 0);

mua{1} = a.MUA(a.DETS(:, 1) == 6, :);
mua{2} = a.MUA(a.DETS(:, 1) == 7, :);

% separate the time series based on whether the stimulus is in the
% receptive field

t_pk  = (a.TD>0) & (a.TD < 0.2);
t_npk = (a.TD > -0.1) & (a.TD < 0);

% reduce the baseline to 0 for each electrode:

elecIdx = [1 : 1210; 1211 : 2420];

idx_in  = {};
idx_out = {};

for k = 1 : 2
    mua{k} = mua{k} - mean(mua{k}(:, t_idx), 2);
    % find response correspond to stimulus inside and outside of the
    % receptive field
    idx_in{k}  = [];
    idx_out{k} = [];
    
    for k1 = 1 : n_locations
        idx = find(a.DETS(elecIdx(1, :), 2) == k1);
        mloc_ts{k}(k1, :) = mean(mua{k}(idx, :));
        if median(mloc_ts{k}(k1, t_pk)) > criterion*abs(median(mloc_ts{k}(k1, t_npk))),
            idx_in{k} = [idx_in{k}, k1];
        else
            idx_out{k} = [idx_out{k}, k1];
        end
        in{k}  = mloc_ts{k}(idx_in{k}, :);
        out{k} = mloc_ts{k}(idx_out{k}, :);
    end
end

% compute mean response
m_in = [];
t_select = (a.TD>-0.1)& (a.TD<=0.3);

for k = 1 : 2
    tmp = mean(in{k});
    tmp = tmp - mean(tmp(t_npk));
    m_in(k, :) = norm_max(tmp(t_select));
end

%%
figure (1), clf

for k = 1 : 2
    subplot(1, 2, k)
    plot(a.TD(a.TD>-0.1), mean(in{k}(:, a.TD>-0.1))), hold on
    plot(a.TD(a.TD>-0.1), mean(out{k}(:, a.TD>-0.1)))
    plot([-0.1, 0.4], [0, 0], 'k:'),
   % ylim([-0.1, 0.2])
end

figure (2), clf
for k = 1 : 2
   subplot(1, 2, k)
   plot(a.TD(t_select), m_in(k, :), 'k-'), hold on
   plot([-0.1, 0.3], [0, 0], 'k:')
end
    

%% fit model to the selected mean response
newt = a.TD(t_select);
stim = zeros(1, length(newt));

stim((newt >= 0)& (newt<= 0.2)) = 1;
t    = [1 : length(stim)]./1000;

lb = [0.07, 0.07, 0.5, 0.001, 0];
ub = [1, 1, 5, 1, 1];

seed = [0.1,  0.1, 2, 0.02, 0.02]; % tau1, tau2, n, sigma, shift, scale

for k = 1 : 2
    param(k, :) = fminsearchbnd(@(x) trf_dcts_finefit(x, m_in(k, :), stim, t), seed, lb, ub);
    ext_param(k, :) = [param(k, 1), 0, param(k, 2), param(k, 3), param(k, 4), param(k, 5), 1];
    pred(k, :) = norm_max(trf_dCTSmodel(ext_param(k, :), stim, t));
end


%% plot stimulus and model fit

figure (2), 
for k = 1 : 2
    subplot(1, 2, k)
    plot(newt, stim, '-', 'color', [.5, .5, .5])
    plot(newt, pred(k, :), 'r-', 'linewidth', 3), box off
end


%% load single unit data
b = xlsread('/Volumes/server/Projects/Temporal_integration/data/ECoG/Figure1_ACE_data.xlsx');

cell{1} = b(2 : 27, 2 : 11);
cell{2} = b(30 : 50, 2 : 11);
cell{3} = b(53 : end, 2 : 11);

% cell1 = b(2 : 27, 2 : 11);
time{1} = b(2 : 27 ,1);
time{2} = b(30 : 50, 1);
time{3} = b(53 : end, 1);

for k = 1 : 3
    rsp{k}  = norm_max(cell{k}(:, 10));
end
stim = [ones(1, 200), zeros(1, 300)];
t    = [1 : length(stim)]./1000;

init{1} = [0.07, 0.07, 2, 0.05, 0.01]; init{2} = init{1};
init{3} = [0.07, 0.1, 2, 0.01, 0.01];

low  = [0, 0, 0, 0, 0];
high = [1, 1, 10, 1, 0.1];
for k = 1 : 3
    s_param(k, :) = fminsearchbnd(@(x) trf_fitModel_contrast(x, rsp{k}', time{k}, 1, stim), init{k}, lb, ub);
    ext_sparam(k, :) = [s_param(k, 1), 0, s_param(k, 2), s_param(k, 3), s_param(k, 4), s_param(k, 5), 1];
    sing_pred{k} = norm_max(trf_dCTSmodel(ext_sparam(k, :), stim, t));
    sing_pred{k} = sing_pred{k}(time{k});
end

%% 
figure (3), clf
for k = 1 : 3
   subplot(3, 1, k)
   plot(time{k}./1000, rsp{k}, 'k.-', 'markersize', 20), hold on
   plot(time{k}./1000, sing_pred{k}, 'r.:', 'markersize', 20), axis tight, box off
   set(gca, 'ytick', [0, 0.5, 1]),
   xlim([0, 1])
end

%% compare the fitted parameters

ecogFile = fullfile(temporalRootPath, 'results', 'ecog1_Output');
b = load(ecogFile);

for k = 1 : 8
    idx = (k - 1) * 100 + 1 : k* 100;
    mfparam(k, :) = b.prm.ecog.model.mParam{k};
end

figure (4), clf
plot(mfparam(:, 3), 'o')

%% compute time to peak and asymptote for each measurement

for k = 1 : 3
    [t2pk(k), asymp(k)] = trf_reparameterize(s_param(k, :));
end
for k = 1 : 2
     [t2pk1(k), asymp1(k)] = trf_reparameterize(param(k, :));
end

%% 

figure (4), clf
subplot(1, 2, 1)
plot(1, mean(t2pk), 'ko', 2, mean(t2pk1), 'ko'), hold on, xlim([0.5, 2.5])
plot([1, 1, 1], t2pk, 'r*'), plot([2, 2], t2pk1, 'r*'), ylim([50, 180])

subplot(1, 2, 2)
plot(1, mean(asymp), 'ko', 2, mean(asymp1), 'ko'), hold on, xlim([0.5, 2.5])
plot([1, 1, 1], asymp, 'r*'), plot([2, 2], asymp1, 'r*'), ylim([0, 0.2])

%% load original prediction
c = load(fullfile(temporalRootPath, 'output', 'df_results_00.mat'));

for k = 1 : 6
   idx = (k - 1) * 100 + 1 : k * 100;
   mparam(k, :) = median(c.rs.ecog1.fprm(idx, :));
   mpred(k, :)  = median(c.rs.ecog1.fpred(idx, :));
end

figure
subplot(1, 3, 1)
plot(mparam(:, 1), 'o')
subplot(1, 3, 2), plot(mparam(:, 2), 'o')
subplot(1, 3, 3), plot(mparam(:, 3), 'o')

%% single unit predict to 500 ms stimulus
stim = [zeros(1, 200), ones(1, 500), zeros(1, 500)];
t = [1 : length(stim)]./1000;

m_sparam = mean(ext_sparam);
m_param  = mean(ext_param);

pred_sing = norm_max(trf_dCTSmodel(m_sparam, stim, t));
pred_mua  = norm_max(trf_dCTSmodel(m_param, stim, t));



figure (5), clf
subplot(2, 1, 1)
plot(t, norm_max(mpred(1, 1 : 1200)), 'k-'), hold on
plot(t, stim, 'k-', t, pred_sing, 'r-'), axis tight
subplot(2, 1, 2)
plot(t, norm_max(mpred(3, 1 : 1200)), 'k-'), hold on
plot(t, stim, 'k-', t, pred_mua, 'r-'), axis tight
