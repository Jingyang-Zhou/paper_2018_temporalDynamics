% simulate Sit et al. 2009

% build a 2-dimensional receptive field

sz = 50;

% make receptive field, in unit of mm
x = linspace(-2.75, 2.75, sz);
y = linspace(-2.75, 2.75, sz);
[X, Y] = meshgrid(x, y);

mu    = [0, 0];
sigma = [1, 0; 0, 1];
rf_sp = mvnpdf([X(:), Y(:)], mu, sigma);
rf_sp = reshape(rf_sp, length(x), length(y));

% add the time dimension to the spatial receptive field
t  = 0 : 0.005 : 0.5;
irf = gammaPDF(t, 0.05, 2);

% make stimulus
stim = zeros(size(rf_sp));
stim(22 : 28, 22 : 28) = 1;

% visualize spatial stimulus and receptive field
figure (1), clf
ax1 = subplot(3, 5, 1); imagesc(stim), axis square, colormap(ax1, 'gray')
set(gca, 'xtick','', 'ytick', ''), title('stimulus contrast'), set(gca, 'fontsize', 12)

ax2 = subplot(3, 5, 3)
imagesc(x, y, rf_sp); title('2d population receptive field'), colormap(ax2, 'jet')
xlabel('distance from the RF center (mm)'), ylabel('distance from the R2 center (mm)')
set(gca, 'fontsize', 12), axis square

subplot(3, 5, 5)
plot(t, irf, 'k-', 'linewidth', 3), box off, set(gca, 'yticklabel', '')
xlabel('time (s)'), title('Impulse response function')

%% add the time dimension to both the stimulus and the impulse

% prep to add time to stimulus
stimLth        = 0.2;
stim_t         = zeros(length(t));

stim_t(t<=stimLth) = 1;
stim_3d        = zeros([size(stim), length(t)]);

% prep to add time to receptive field
rf = zeros([size(rf_sp), length(t)]);

% add the time dimension to both stimulus and receptive field
for k = 1 : length(t)
   stim_3d(:, :, k) = stim.*stim_t(k); 
   rf(:, :, k)      = rf_sp.*irf(k);
end

%% static nonlinearity models

% where the model fails according to Sit et al. 2009:
% (1) spatial profile of the response shouldn't widen and change shape as
%     stimulus contrast increases
% (2) The model shouldn't predict a longer latency for the falling edge of
%     the response at high contrast

% make population prediction to a low and a hight contrast stimulus
ctr_levels = [0.2, 1];
stim_ctr{1} = stim_3d.*ctr_levels(1);
stim_ctr{2} = stim_3d.*ctr_levels(2);

% Step 1 : linear computation
rsp_static = {};

for k = 1 : 2
    % static nonlinear response
   tmp           = convn(stim_ctr{k}, rf, 'full'); 
   rsp_linear{k} = tmp(26:75, 26:75, 1 : length(t));
   rsp_static{k} = rsp_linear{k}.^0.3;
end

% normalize the responses
rsp_max = max(rsp_static{2}(:));
for k = 1 : 2
    rsp_static{k} = rsp_static{k}./rsp_max;
end

%% static nonlinearity models: visualization

figure (1),
idx = [1, 3]; title_txt = {'low contrast spatial response', 'high contrast spatilal response'};
for k = 1 : 2
    subplot(3, 5, 5 + idx(k)),colormap('jet')
    imagesc(x, y, squeeze(rsp_static{k}(:, :, 40)), [0, 1.2]), hold on, axis square
    % outline the equi-contour
    contour(X, Y, squeeze(rsp_static{k}(:, :, 40)), 'w')
    plot(x(25)*ones(1, length(y)), y, 'k-', 'linewidth', 3), title(title_txt{k}), set(gca, 'fontsize', 12)
end

% take a spatial slice and visualize the time component
for k = 1 : 2
    subplot(3, 5, 5 + idx(k)+1),cla, colormap('jet'), set(gca, 'colorOrder', parula(sz))
    patch([0, stimLth, stimLth, 0], [0, 0, ctr_levels(k), ctr_levels(k)], 0.5 * ones(1, 3)), hold on
    plot(t, squeeze(rsp_static{k}(25, :, :))'), box off, axis tight, xlabel('time (s)'), 
    set(gca, 'ytick', ''), ylim([0, max(rsp_static{2}(:))]), set(gca, 'fontsize', 12)
end

% compare one time course for low and for high contrast stimulus
subplot(3, 5, 10), cla, set(gca, 'colororder', gray(3))
for k = 1 : 2
    rsp_tmp = squeeze(rsp_static{k}(25, 25, :));
    rsp_tmp = rsp_tmp / max(rsp_tmp);
    plot(t, rsp_tmp, '--', 'linewidth', 3), hold on
    %plot(t, squeeze(rsp_static{k}(25, 25, :)), '-', 'linewidth', 3), hold on
    set(gca, 'fontsize', 12)
end
axis tight, box off, xlabel('time (s)'), legend('rsp to low ctr', 'rsp to high ctr')

%% lateral propogation models

v = 1; % conduction speed
c = [x(round(sz/2)), y(round(sz/2))]; % center of the kernal

%% normalization models

% linear sum of the response pool was computed before, as rsp_linear
% linear pool receptive field:
% mu    = [0, 0];
% sigma = [1, -0.15; -0.15, 1];

% First compute the linear sum from the normalization pool:
% receptive field for the normalization pool:
mu_norm    = [0, 0];
sigma_norm = [1.5, -0.15; -0.15, 1.5];

% build normalization receptive field
rf_n_sp = mvnpdf([X(:), Y(:)], mu, sigma);
rf_n_sp = reshape(rf_sp, length(x), length(y));

% add temporal element to the receptive field:
