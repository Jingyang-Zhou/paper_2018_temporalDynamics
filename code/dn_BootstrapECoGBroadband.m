% dn_bootstrap ecog data

function bs_bbts = dn_BootstrapECoGBroadband(bbts, t, nBoots)
%% FIGURE KNOB

figureOn = 0;

%% EXAMPLES

% dataLoc    = fullfile(dn_ECoG_RootPath, 'data');
% dt_fName   = 'dn_preprocessedData.mat';
% prm_fName  = 'dn_params.mat';
%
% a = load(fullfile(dataLoc, dt_fName));
% b = load(fullfile(dataLoc, prm_fName));
%
%  bbts   = a.dt.ecog.bbts_roi;
%  t      = b.prm.ecog.t;
%  stim   = b.prm.ecog.stim;
%
%  nBoots = 100;

%% DERIVED PARAMETERS

nrois = length(bbts);
t_lth = length(t);

%% BOOTSTRAP THE MEAN OF THE DATA FROM EACH ROI

for iroi = 1 : nrois
    nElecPerRoi   = size(bbts{iroi}, 2);
    if nElecPerRoi ~= 0,
        for iBoot = 1 : nBoots
            idx  = randi(nElecPerRoi, [1, nElecPerRoi]); % bootstrap sample selection
            count = (iroi - 1) * nBoots + iBoot;
            bs_bbts(:, count) = mean(bbts{iroi}(:, idx), 2);
        end
    else
        bs_bbts(:, (iroi - 1) * nBoots + 1 : iroi * nBoots) = nan;
    end
end

%% VISUALIZE

if figureOn
    % USEFUL FUNCTIONS ---------------------------------------------------
    normMax = @(x) x./max(x);
    % PLOTTING -----------------------------------------------------------
    figure, % plot the roi series without normalization to the max
    for iroi = 1 : nrois
        subplot(2, 3, iroi),
        count = (iroi - 1) * nBoots + 1 : iroi * nBoots;
        
        % plot mean
        m = mean(bs_bbts(:, count), 2);
        s = std(bs_bbts(:, count), [], 2);
        shadedErrorBar(t, m, s, 'k-'),
        xlim([0, max(t)]), ylim([-1, 8])
    end
    
    figure, set(gca, 'colorOrder', copper(6)), hold on,  % plot the roi series with normalization
    for iroi = 1 : nrois
        count = (iroi - 1) * nBoots + 1 : iroi * nBoots;
        m = normMax(mean(bs_bbts(:, count), 2));
        plot(t, m, 'linewidth', 2), xlim([0, max(t)])
    end
end