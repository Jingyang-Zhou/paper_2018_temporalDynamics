%% Simulate Sit et al. 2009

%% Some thoughts about the simulation

% (1) In the JNS paper, instead of a naka rushton function, we assumed a
% simpler shape of the CRF. Assumption 1 is that we bypass the small
% contrast range in the fMRI experiment.

% (2) In the paper, the response time course is predicted by concatenating 2
% sigmoidal functions(convoled with stimulus time course?) with different lag parameters. 

%%
sz        = 50; % number of entries per dimension/how dense the samples are
sp_width  = 5; % in unit of mm
sigma_val = 2; % receptive field size, or the size of the smallest receptive field

%% make 2-d receptive fields

% RF center, assume equally spaced RF centers
x      = linspace(-sp_width, sp_width, sz); 
y      = x; 
[X, Y] = meshgrid(x, y);

% parameters for the linear pool
sigma_lin = [1, 0; 0, 1];
mu    = [0, 0];
rf_lin_sp = mvnpdf([X(:), Y(:)], mu, sigma_lin);
rf_lin_sp = reshape(rf_lin_sp, length(x), length(y));

% parameters for the normalization pool
sigma_norm = [1.5, 0; 0, 1.5];
rf_norm_sp = mvnpdf([X(:), Y(:)], mu, sigma_norm);
rf_norm_sp = reshape(rf_norm_sp, length(x), length(y));

% contrast response function
contrast = 0 : 0.1 : 1;
prm.dm = 2; prm.s50 = 0.3; prm.n = 4;
crf    = nakaRushton(prm, contrast);

%% make stimulus

% spatial aspect of the stimulus
stim_s = zeros(size(rf_lin_sp));

center  = round(sz/2);
stim_sz = 2; % pixels
stim_s(center-stim_sz : center+stim_sz, center-stim_sz : center+stim_sz) = 1;

% stimulus time course and the 3-d stimulus
t  = 0 : 0.005 : 0.5; t_lth = length(t);
stim_t = zeros(1, length(t));
stim_t((t > 0.1) &(t < 0.3)) = 1;

% make 3-d stimulus
Ix = zeros([size(stim_s), t_lth]);
Ix(:, :, stim_t == 1) = repmat(stim_s, [1, 1, nnz(stim_t)]);

%% static nonlinear model

% impulse repsonse function in time
irf = gammaPDF(t, 0.05, 2);

%% lateral propogation model

%% normalization model - onset
stim_case = 'onset'; % options: 'onset', 'steadystate', 'offset' 
onset_idx = 1 : 50;

low_contrast = 0.5;

% model parameters:
g0 = 0.3;
C  = 0.5;

% high contrast
Vx.onset = sit2009_normalization(Ix(:, :, onset_idx), t(onset_idx), stim_t(onset_idx), ...,
    g0, C, rf_lin_sp, rf_norm_sp, stim_case);
% low contrast
Vx.lowctr_onset = sit2009_normalization(Ix(:, :, onset_idx).*low_contrast, t(onset_idx), stim_t(onset_idx), ...,
    g0, C, rf_lin_sp, rf_norm_sp, stim_case);

%% normalization model - steady state

stim_case = 'steadystate'; % options: 'onset', 'steadystate', 'offset' 
ss_idx    = 22:60;

Vx.ss = sit2009_normalization(Ix(:, :, ss_idx), t(ss_idx), stim_t(ss_idx), ...,
    g0, C, rf_lin_sp, rf_norm_sp, stim_case);

Vx.lowctr_ss = sit2009_normalization(Ix(:, :, ss_idx).*low_contrast, t(ss_idx), stim_t(ss_idx), ...,
    g0, C, rf_lin_sp, rf_norm_sp, stim_case);

%% normalization model - offset

stim_case  = 'offset'; % options: 'onset', 'steadystate', 'offset' 
offset_idx = 40 : 100;

Vx.offset = sit2009_normalization(Ix(:, :, offset_idx), t(offset_idx), stim_t(offset_idx), ...,
    g0, C, rf_lin_sp, rf_norm_sp, stim_case);

Vx.lowctr_offset = sit2009_normalization(Ix(:, :, offset_idx).*low_contrast, t(offset_idx), stim_t(offset_idx), ...,
    g0, C, rf_lin_sp, rf_norm_sp, stim_case);

%% Visualize normalization responses

figure (1), clf
subplot(3, 5, 1),  % plot the relevant stimulus time course : onset
plot(t, stim_t, 'k-'), hold on, plot(t(onset_idx), stim_t(onset_idx), 'r-', 'linewidth', 3), 
box off, xlabel('time (s)'), ylabel('contrast')

subplot(3, 5, 6),  % plot the relevant stimulus time course : steady state
plot(t, stim_t, 'k-'), hold on, plot(t(ss_idx), stim_t(ss_idx), 'r-', 'linewidth', 3)
box off, xlabel('time (s)'), ylabel('contrast')

subplot(3, 5, 11), % plot the relevant stimulus time course : offset
plot(t, stim_t, 'k-'), hold on, plot(t(offset_idx), stim_t(offset_idx), 'r-', 'linewidth', 3)
box off, xlabel('time (s)'), ylabel('contrast')



% difference between high and low contrast response
subplot(3, 5, 2), imagesc(x, y, squeeze(Vx.onset(:, :, 25)), [0 .3]), hold on, axis square,
contour(X, Y, squeeze(Vx.onset(:, :, 25)), 'w-')

subplot(3, 5, 3), imagesc(x, y, squeeze(Vx.lowctr_onset(:, :, 25)), [0 .3]), hold on,  axis square
contour(X, Y, squeeze(Vx.lowctr_onset(:, :, 25)), 'w-')

subplot(3, 5, 7), imagesc(x, y, squeeze(Vx.ss(:, :, 25)), [0, 0.3]), hold on, axis square,
contour(X, Y, squeeze(Vx.ss(:, :, 25)), 'w-')

subplot(3, 5, 8), imagesc(x, y, squeeze(Vx.lowctr_ss(:, :, 25)), [0, 0.3]), hold on,  axis square
contour(X, Y, squeeze(Vx.lowctr_ss(:, :, 25)), 'w-')

subplot(3, 5, 12), imagesc(x, y, squeeze(Vx.offset(:, :, 1))), hold on, axis square,
contour(X, Y, squeeze(Vx.offset(:, :, 1)), 'w-')

subplot(3, 5, 13), imagesc(x, y, squeeze(Vx.lowctr_offset(:, :, 1))), hold on,  axis square
contour(X, Y, squeeze(Vx.lowctr_offset(:, :, 1)), 'w-')



% plot a spatial slice over time - onset
subplot(3, 5, 4), set(gca, 'ColorOrder', parula(length(onset_idx))), hold on
plot(t(onset_idx), squeeze(Vx.onset(25, :, :))'), box off

subplot(3, 5, 5), set(gca, 'ColorOrder', parula(length(onset_idx))), hold on
plot(x, squeeze(Vx.onset(25, :, :))), box off

% plot a spatial slice over time
subplot(3, 5, 9), set(gca, 'ColorOrder', parula(length(onset_idx))), hold on
plot(t(onset_idx), squeeze(Vx.ss(25, :, :))'), box off

subplot(3, 5, 10), set(gca, 'ColorOrder', parula(length(onset_idx))), hold on
plot(x, squeeze(Vx.ss(25, :, :))), box off


% plot a spatial slice over time
subplot(3, 5, 14), set(gca, 'ColorOrder', parula(length(offset_idx)))
plot(t(onset_idx), squeeze(Vx.offset(25, :, :))'), box off

subplot(3, 5, 15), set(gca, 'ColorOrder', parula(length(offset_idx)))
plot(x, squeeze(Vx.offset(25, :, :))), box off

