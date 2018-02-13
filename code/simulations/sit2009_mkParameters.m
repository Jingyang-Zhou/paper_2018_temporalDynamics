% simulate Sit et al. 2009: make parameter files

function prm = sit2009_mkParameters(sz, rf_sig, tau, cort_sz, stim_dur)

%% Test

figureOn = 0;

% sz       = 30;
% rf_sig   = 1;
% tau      = 0.05;
% cort_sz  = 2.75; % mm, the extent of cortex we are looking at
% stim_dur = 0.2;

%% (1.1)initiate parameters that apply to both the cortex and the stimulus

prm = [];

prm.t = 0 : 0.05 : 1; t_lth = length(prm.t);

%% (2.1) MAKE RECEPTIVE FIELD : INITIATE
prm.rf = [];

% number of grid point/receptive field centers on the cortex
prm.rf.sz = sz;

% receptive field parameters: center and spread
prm.rf.mu    = [0, 0];
prm.rf.sigma = [rf_sig, 0; 0, rf_sig];

%% (2.2) MAKE RECEPTIVE FIELD : SPATIAL

prm.rf.x = linspace(-cort_sz, cort_sz, sz);
prm.rf.y = prm.rf.x; % assuming the ROI is a square in millimeters

[prm.rf.X, prm.rf.Y] = meshgrid(prm.rf.x, prm.rf.y);

prm.rf.f_s = mvnpdf([prm.rf.X(:), prm.rf.Y(:)], prm.rf.mu, prm.rf.sigma);
prm.rf.f_s = reshape(prm.rf.f_s, sz, sz);

%% (2.3) MAKE RECEPTIVE FIELD : TEMPORAL

prm.rf.f_t = gammaPDF(prm.t, tau, 2);

%% (2.4) MAKE RECEPTIVE FIELD : 3-D

prm.rf.f = zeros([size(prm.rf.f_s), length(prm.rf.f_t)]);

for it = 1 : t_lth
    prm.rf.f(:, :, it) = prm.rf.f_s.*prm.rf.f_t(it);
end

%% (3.1) MAKE STIMULUS : THE SPATIAL EXTENT

prm.stim = [];

prm.stim.s = zeros(size(prm.rf.f_s));

center = round(sz)/2;
spread = 3; % number of pixels

prm.stim.s(center - spread : center + spread, center - spread : center + spread) = 1;

%% (3.2) MAKE STIMULUS : TIME COURSE

prm.stim.t = zeros(1, t_lth);
prm.stim.t(prm.t > 0.1 & prm.t < stim_dur + 0.1) = 1;

%% (3.3) MAKE STIMULUS : 3D

prm.stim.st = zeros([size(prm.stim.s), length(prm.stim.t)]);

prm.stim.st(:, :, prm.stim.t == 1) = repmat(prm.stim.s, [1, 1, nnz(prm.stim.t)]);

%% (4) VISUALIZE

if figureOn
    figure (100), clf
    for it = 1 : t_lth
        subplot(1, 2, 1), cla, colormap gray
        imagesc(prm.rf.x, prm.rf.y, prm.rf.f(:, :, it), [0, 0.1]),
        subplot(1, 2, 2), cla, colormap gray
        imagesc(prm.stim.st(:, :, it), [0, 1]), pause(0.2)
        drawnow
    end
end

%%

end