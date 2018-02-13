function rsp = sit2009_staticNonlinear(sz, stim, rf, modelType)

% QUESTION: 

% About modelType:
% (1) "STN" - supressive temporal normalization?
% (2) "ETC" - exponentiated temporal compressive?

%% test:

figureOn = 0;

% sz       = 30;   % number of entries per dimension of receptive field and stimulus
% rf_sig   = 1;    % the spread of the spatial receptive fielf
% tau      = 0.05; % the length of the temporal receptive window
% cort_sz  = 2.75; % mm, the extent of cortex we are looking at
% stim_dur = 0.2;  % stimulus duration
% 
% % make parameters
% prm = sit2009_mkParameters(sz, rf_sig, tau, cort_sz, stim_dur);
% 
% rf   = prm.rf;
% stim = prm.stim;
% 
% modelType = 'STN';

%% construct response output

% useful functions
normMax = @(x) x./max(x(:));

% compute STN model prediction
stn = @(x) x.^2./(0.1^2 + x.^2);

% compute ETC model prediction
etc = @(x) x.^0.4;

rsp = [];

% compute linear response
rsp.lin = convn(stim.st, rf.f, 'full');
mid     = round(size(rsp.lin, 1)/2);
rsp.lin = rsp.lin(mid-sz/2 : mid + sz/2-1, mid-sz/2 : mid + sz/2-1, 1 : size(rf.f, 3));

switch modelType
    case 'STN'
        % compute STN response
        rsp.s_nlin = stn(rsp.lin);
    case 'ETC'
        % compute ETC response
        rsp.s_nlin = etc(rsp.lin);
end

%% visualize

if figureOn
    l_max = max(rsp.lin(:));
    
    for it = 1 : size(rsp.lin, 3)
        figure (1), clf
        subplot(1, 2, 1), imagesc(rsp.lin(:, :, it), [0, l_max])
        subplot(1, 2, 2), imagesc(rsp.s_nlin(:, :, it), [0, 1]), pause(0.5)
    end
end

end