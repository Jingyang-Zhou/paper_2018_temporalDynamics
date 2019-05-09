t = 0.001 : 0.001 : 1;




stim = zeros(5, length(t));

onset = 201;
dur   = 60;
stim(1, onset) = 1;                             % impulse
stim(2, onset+(1:dur)) = 1;                     % calibration (don't plot!)
stim(3, onset+(1:dur*8)) = 1;                   % long pulse
stim(4, [onset*3+(1:dur*2)]) = 1;               % twice the duration
stim(5, [onset+(1:dur) onset*2+(1:dur)]) = 1;   % two pulses
stim(6, :) = stim(3, :) * .2;                  % low contrast
%stim(6, onset+(1:dur*3)) = .2;                  % low contrast

stim = stim(:,1:length(t));

str{1} = 'Impulse Response Function';
str{2} = 'Calibration. DON''T PLOT!';
str{3} = 'Response Reduction for Prolonged Stimulus';
str{4} = 'Sub-additive Temporal Sumamtion';
str{5} = 'Response Reduction for Repeated Stimulus';
str{6} = 'Different Dynamics at Low Contrast';
% 
% xt(1).t = onset;
% xt(1).label = 'Impulse';
% 
% xt(3).t     = [onset+dur/2 onset*3+dur];
% xt(3).label = {'1x' '2x'};
% 
% xt(4).t     = [onset+dur/2 onset*3+dur];
% xt(4).label = {'Stim 1' 'Repated Stim'};


% DN predictions
param = [0.07, 0, 0.1, 2, 0.1, 0, 1];
param = [0.07, 0, 0.1, 2, 0.1, 0, 1];
%param = [0.07, 0, 0.15, 2.2, 0.1, 0, 1];
%param = [0.08, 0, 0.18, 2.2, 0.1, 0, 1];

predDn = dn_DNmodel(param, stim, t);

% normalize
predDn = predDn / max(predDn(2,:));

% cal_stim = stim(2,:);
% cal_dn = predDn(2,:);
% 
% % fit 
% tau = .009; n = 5; delay = .010;
% 
% x0 = [tau n delay];
% lb = [.001 1 0];
% ub = [.1 10 .015];
% 
% 
% irfFun = @(x) fitIRFtoDN(t,cal_stim, x)-cal_dn;
% 
% %irfFun = @(x) convCut(cal_stim, gammaPDF_JW(t, x), length(t))-cal_dn;
% options = optimoptions(@lsqnonlin,'Display','iter');
% 
% params = lsqnonlin(irfFun,x0, lb, ub, options);
% 
% irfFit = gammaPDF_JW(t, params);
% 
% lin = convCut(cal_stim, irfFit, length(t));
% lin = lin / max(lin);
% 
% figure(1), clf; plot(t, cal_dn, t,lin);


irf1   = gammaPDF(t, 0.029, 2); dt = 14;
%irf1   = gammaPDF(t, 0.024, 2); dt = 18;
irf1   = [zeros(1,dt) irf1(1:end-dt)];

lin1 = conv2(stim, irf1, 'full');
lin1 = lin1(:, 1:length(t));


lin1 = lin1 / max(lin1(2,:));




fH = figure(2); clf; set(fH, 'Color', 'w', 'Position', [1 1 1000 1200])

% Impulse response function in first column
subplot(3,2,3)
set(gca, 'FontSize', 16, 'XTick', [], 'YTick', [], 'LineWidth', 2); 
   hold on;box off;
plot(t, stim(1,:)/50, 'k-',  t, lin1(1,:), 'k:', 'LineWidth', 2);
title(str{1})

% Non linear phenomena in second column
for ii = 3:6
   subplot(4,2,2*(ii-2));
   set(gca, 'FontSize', 15, 'XTick', [],  'LineWidth', 2, 'YTick', [])
   hold on;
   
   if ii == 4
       thisstim =  stim(ii,:)   + stim(2,:);
       thisPred =  predDn(ii,:) + predDn(2,:);
       thisLin  =  lin1(ii,:)   + lin1(2,:);
   else
       thisstim =  stim(ii,:);
       thisPred =  predDn(ii,:);
       thisLin  =  lin1(ii,:);
   end
   
   fill(t, thisstim, .7*[1 1 1])
   plot(t, thisstim, 'k-', t, thisPred, 'r-', t, thisLin, 'k:', 'LineWidth', 2);%, t, lin2(ii,:), 'g--')
   
   %
   % fill(t, stim(ii,:), .7*[1 1 1])
   % plot(t, stim(ii,:), 'k-', t, predDn(ii,:), 'r-', t, lin1(ii,:), 'k:', 'LineWidth', 2);%, t, lin2(ii,:), 'g--')
   
   
   ylim([0 1.7])
   box off;
   title(str{ii})
   
end
    
hgexport(fH, '~/Desktop/temporalPhenomena.eps')