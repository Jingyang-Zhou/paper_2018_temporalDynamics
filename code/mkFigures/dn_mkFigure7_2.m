% compare DN and 2 temporal channels model

saveFigure = 0;

%% PRE-DEFINED PARAMETERS

% make DN model:
dn_param = [0.08, 0, 0.08, 2, 0.1, 0, 1];

% make 2TC mode:
ttc_param.b1 = 0.3; ttc_param.b2 = 0.1;

% time course
dt = 0.001;
t  = dt : dt : 1.2;

%% PREDICTION TO A HIGH AND A LOW CONTRAST STIMULUS


% MAKE STIMULUS -----------------------------------------------------------
stim = [];

stim.highCtrst = zeros(1, length(t));
stim.highCtrst(t>0.2 & t <= 0.7) = 1;

lowCtrst      = 0.2;
stim.lowCtrst = stim.highCtrst .*lowCtrst;

% figure (1), clf
% subplot(1, 2, 1), plot(t, stim.highCtrst, 'k-'), hold on, ylim([0, 1]), box off
% subplot(1, 2, 2), plot(t, stim.lowCtrst, 'k-'), hold on, ylim([0, 1]), box off

% MAKE MODEL PREDICTIONS --------------------------------------------------

% DN MODEL:
dnPred(1, :)  = dn_DNmodel(dn_param, stim.highCtrst, t);
dnPred(2, :)  = dn_DNmodel(dn_param, stim.lowCtrst, t);
dnPred        = dnPred ./ max(dnPred(:));

% 2TC MODEL:
ttcPred(1, :) = dn_2Chansmodel(ttc_param, stim.highCtrst, t, dt);
ttcPred(2, :) = dn_2Chansmodel(ttc_param, stim.lowCtrst, t, dt);
ttcPred       = ttcPred ./ max(ttcPred(:));

figure (1), clf
for k = 1 : 2
    subplot(1, 3, k)
    plot(t-0.2, stim.highCtrst, 'k-', 'linewidth', 3), hold on
    plot(t-0.2, stim.lowCtrst, '-', 'color', (1-lowCtrst) * ones(1, 3), 'linewidth', 3),
    box off
end
subplot(1, 3, 1), plot(t-0.2, dnPred', 'r-', 'linewidth', 3), ylim([-0.1, 1])
index = find(dnPred(1, :) == max(dnPred(1, :))); new_t = t - 0.2;
plot([new_t(index), new_t(index)], [0, 1], 'r:', 'linewidth', 2)

subplot(1, 3, 2), plot(t-0.2, ttcPred', 'b-', 'linewidth', 3), ylim([-0.1, 1])
index = find(ttcPred(1, :) == max(ttcPred(1, :))); new_t = t - 0.2;
plot([new_t(index), new_t(index)], [0, 1], 'b:', 'linewidth', 2)

for k = 1 : 2
    subplot(1, 3, k)
    set(gca, 'xtick', [0, 0.5], 'ytick', [0, lowCtrst, 1], 'fontsize', 18)
    xlabel('time (s)'), ylabel('amplitude'), axis tight
end

%% CARANDINI ET AL. 1997

c97_param = [0.05, 1]; % tau1 and w
rcprm.C  = 5; % capacitance
rcprm.g0 = 0.1;
rcprm.k  = 1; % determines the effectiveness of the normalization pool

c97pred(1, :) = dn_Carandini1997(c97_param, rcprm, t, stim.highCtrst);
c97pred(2, :) = dn_Carandini1997(c97_param, rcprm, t, stim.lowCtrst);
c97pred       = abs(c97pred./max(c97pred(:)));

figure (1), 
subplot(1, 3, 3)
plot(t-0.2, stim.highCtrst, 'k-', 'linewidth', 3), hold on, box off
plot(t-0.2, stim.lowCtrst, '-', 'color', (1-lowCtrst) * ones(1, 3), 'linewidth', 3),
plot(t-0.2, c97pred, 'r-', 'linewidth', 3)

index = find(c97pred(1, :) == max(c97pred(1, :))); new_t = t - 0.2;
index1 = find(c97pred(2, :) == max(c97pred(2, :))); new_t = t - 0.2;
plot([new_t(index), new_t(index)], [0, 1], 'b:', 'linewidth', 2)
plot([new_t(index1), new_t(index1)], [0, max(c97pred(2, :))], 'b:', 'linewidth', 2)

set(gca, 'xtick', [0, 0.5], 'ytick', [0, lowCtrst, 1], 'fontsize', 18)
xlabel('time (s)'), ylabel('amplitude'), axis tight

%% SAVE FIGURE

if saveFigure
    figureLoc = fullfile(dn_ECoG_RootPath, 'analysisFigures');
    printnice(1, 0, figureLoc, 'fg_compareDN_2TC_tmp')
    %printnice(2, 0, figureLoc, 'fg_compareDN_c97')
end



