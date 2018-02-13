function rsp = sit2009_DN(sz, stim, rf, t, modelType)

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
% rf    = prm.rf;
% stim  = prm.stim;
% t_lth = length(prm.t);
% 
% modelType = 'monophasic'; % monophasic or biphasic
% 
%% Initiate: make model parameters

dn.prm = [];

switch modelType
    case 'monophasic'
        dn.prm = [0.07, 0, 0.1, 2, 0.05, 0, 1]; % tau1, weight, tau2, n, sigma, shift, scale
    case 'biphasic'
        dn.prm = [0.05, 0.3, 0.05, 2, 0.05, 0, 1];
end

t_lth = length(t);

%% compute linear response in space

rsp_s = [];
rsp_s = conv2(stim.s, rf.f_s, 'same');

%% compute DN response
rsp = [];

for it = 1 : t_lth
   rsp(:, :, it) = rsp_s.* stim.t(it);
end

rsp = reshape(rsp, [sz^2, t_lth]);
rsp = dn_DNmodel(dn.prm, rsp, t);
rsp = reshape(rsp, [sz, sz, t_lth]);

%% visualize

if figureOn
    l_max = max(rsp(:));
    
    for it = 1 : size(rsp, 3)
        figure (100), 
        imagesc(rsp(:, :, it), [0, l_max]), pause(0.2), drawnow
    end
end

end