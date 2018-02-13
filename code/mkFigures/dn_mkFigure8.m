% dn_mkFigure 8

%% load ECoG data

dataLoc = fullfile(temporalRootPath, 'data');
fName   = 'dn_data.mat';

a = load(fullfile(dataLoc, fName));
goodChans = a.raw.goodChannels;

%% useful functions

normMax_range = @(x, range) (x - mean(x(range)))./max(x - mean(x(range)));
normMax       = @(x) x./max(x);

%% find a V1 and V3 foveal and peripheral electrode

electrodeNum = [92, 68, 115, 104];
% 92 : V2 fovea
% 68 : V2 periphery

% 115 : V3AB fovea
% 104: VO periphery

nElec = length(electrodeNum);

for k = 1 : nElec
   idx = find(goodChans == electrodeNum(k));
   ts(k, :) = normMax_range(mean(a.raw.bbmatrix(:, :, idx), 2), 1 : 200);
end

t = [1 : length(ts)]./1000;

% visualize the time courses

lineType = {'k-.', 'k:'};

figure (1), clf
for k1 = 1 : 3
    subplot_tight(2, 3, k1, 0.02), count = 0;
    for k = [1 3], count = count + 1; plot(t, ts(k, :)', lineType{count},'linewidth', 2), hold on,  end
    patch([0.2, 0.7, 0.7, 0.2], [0, 0, 1, 1], 'k', 'facealpha', 0.2),
    set(gca, 'xaxislocation', 'origin', 'xtick', [0.2, 0.7], 'ytick', [0, 1], 'fontsize', 12)
    axis tight, axis square, box off
end

for k1 = 1 : 3
    subplot_tight(2, 3, k1 + 3, 0.02), count = 0;
    for k = [2 4], count = count + 1; plot(t, ts(k, :)', lineType{count}, 'linewidth', 2), hold on,  end
    patch([0.2, 0.7, 0.7, 0.2], [0, 0, 1, 1], 'k', 'facealpha', 0.2),
    set(gca, 'xaxislocation', 'origin', 'xtick', [0.2, 0.7], 'ytick', [0, 1], 'fontsize', 12)
    axis tight, axis square, box off
end

%% load ECoG data

dataLoc = fullfile(temporalRootPath, 'data');
fName   = 'dn_data.mat';

a       = load(fullfile(dataLoc, fName));
param   = a.param;

%% fit uniphasic and biphasic DN model

stim = zeros(1, length(t)); stim(201 : 700) = 1;

seed1 = [0.1, 0.1, 2, 0.05, 0];
seed2 = [0.05, 0.95, 0.1, 0.05, 0];

% fit uniphasic model
for k = 1 : 4
    [modelprm(k, :), dnprd_uni(k, :), ~, ~] = dn_fineFit(ts(k, :), stim, t, param, seed1, 'uniphasic');
end
% fit biphasic model
for k = 1 : 4
    
    if k == 2, seed2 = [0.05, 0.98, 0.1, 0.1, 0];
    elseif k == 4, seed2 = [0.05, 0.8, 0.1, 0.1, 0];
    else seed2 = [0.05, 0, 0.1, 0.1, 0];
    end
    [modelprm(k, :), dnprd_bi(k, :), ~, ~] = dn_fineFit(ts(k, :), stim, t, param, seed2, 'biphasic');
end

% fit 2tc channel
fields = {'b1', 'b2', 'epsilon'};
dt     = 0.001;
prm = struct();

for k = 1 : 4
    ttcprm(k, :) = fminsearch(@(x) dn_fit2ChansModel(x, ts(k, :), t, stim, dt, fields), [0.3, 0.8, 0.0002]);
    prm.b1 = ttcprm(k, 1);
    prm.b2 = ttcprm(k, 2);
    prm.epsilon = ttcprm(k, 3);
    ttcprd(k, :) = dn_2Chansmodel(prm, stim, t, dt);
end


lineCol = {'r', 'b'};

figure (1), subplot_tight(2, 3, 1, 0.02), count = 0;
% plot dn fit
for k = [1, 3], count = count + 1; plot(t, normMax_range(dnprd_uni(k, :), [1 : 200]), 'linewidth', 2, 'color', lineCol{count}), end
legend('V2 fovea', 'VO fovea', 'stim', 'V2 DN fit', 'VO DN fit')

subplot_tight(2, 3, 4, 0.02), count = 0;
% plot dn fit
for k = [2, 4], count = count + 1; plot(t, normMax_range(dnprd_uni(k, :), [1 : 200]), 'linewidth', 2, 'color', lineCol{count}), end
legend('V2 peri.', 'V3AB peri.', 'stim', 'V2 2TC fit', 'V3ab 2TC fit')

subplot_tight(2, 3, 2, 0.02), count = 0;
% plot 2tc fit
for k = [1, 3], count = count + 1; plot(t, normMax_range(ttcprd(k, :), [1 : 200]), 'linewidth', 2, 'color', lineCol{count}), end
legend('V2 fovea', 'VO fovea', 'stim', 'V2 DN fit', 'VO DN fit')

subplot_tight(2, 3, 5, 0.02), count = 0;
% plot 2tc fit
for k = [2, 4], count = count + 1; plot(t, normMax_range(ttcprd(k, :), [1 : 200]), 'linewidth', 2, 'color', lineCol{count}), end
legend('V2 fovea', 'VO fovea', 'stim', 'V2 DN fit', 'VO DN fit')

subplot_tight(2, 3, 3, 0.02), count = 0;
% plot bi_dn fit
for k = [1, 3], count = count + 1; plot(t, normMax_range(dnprd_bi(k, :), [1 : 200]), 'linewidth', 2, 'color', lineCol{count}), end
legend('V2 fovea', 'VO fovea', 'stim', 'V2 DN fit', 'VO DN fit')

subplot_tight(2, 3, 6, 0.02), count = 0;
% plot bi_dn fit
for k = [2, 4], count = count + 1; plot(t, normMax_range(dnprd_bi(k, :), [1 : 200]), 'linewidth', 2, 'color', lineCol{count}), end
legend('V2 fovea', 'VO fovea', 'stim', 'V2 DN fit', 'VO DN fit')


%% fit 2 temporal channels model

