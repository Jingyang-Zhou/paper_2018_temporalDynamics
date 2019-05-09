% review_contrast simultation

%% MAKE A STIMULUS

t    = 0.001 : 0.001 : .2;
stim = zeros(size(t));
stim(1 : 200) = 1;

% make contrast stimulus

nContrast = 10;
contrast  = linspace(0.1, 1, nContrast);
cStim     = zeros(nContrast, length(t));
for k = 1 : nContrast
   cStim(k, :) = stim .* contrast(k);
end

%% MAKE DN MODEL PREDICTIONS

prd = [];

param = [0.02, 0, 0.03, 2, 0.2, 0, 1];

for k = 1 : nContrast
   prd(k, :) = dn_DNmodel(param, cStim(k, :), t);
end

figure (1), clf
subplot(121), set (gca, 'colororder', copper(nContrast)), hold on, plot(t, cStim, 'linewidth', 2)
subplot(122), set (gca, 'colororder', copper(nContrast)), hold on, plot(t, prd, 'linewidth', 2), 
xlabel('time (s)'), set(gca, 'ytick', '', 'fontsize', 14),

%% GENERATE A NAKA-RUSHTON EQUATION

prm = [];

prm(1, :) = [0.0395, 0, 0.0397, 2.2829, 0.0266, 0, 1];
prm(2, :) = [0.0436, 0, 0.1049, 1.6136, 0.0027, 0, 1];
prm(3, :) = [0.0337, 0, 0.1373, 3.5056, 0.016,  0, 1];

cprd = [];

for k = 1 : nContrast
   for icell = 1 : 3
      cprd(icell, k, :) = dn_DNmodel(prm(icell, :), cStim(k, :), t); 
   end
end

figure (2), clf
for k = 1 : 3
   subplot(2, 3, k), set(gca, 'colororder', copper(nContrast)), hold on
   cprd(k, :, :) = cprd(k, :, :) ./max(max(cprd(k, :)));
   plot(t, squeeze(cprd(k, :, :)), 'linewidth', 2)
end

%% FIT THE SHIFTED AND SCALED COPY MODEL

init = [1, 0.04];
lb   = [0, 0];
ub   = [2, 0.15];

% compute average shape
for k = 1 : 3
    tmp = squeeze(cprd(k,  :, :));
    for k1 = 1 : nContrast
       tmp(k1, :) = tmp(k1, :) ./ max(tmp(k1, :));
    end
   aveshape(k, :) =  mean(tmp);
end

for k = 1 : 3
    example = squeeze(cprd(k, 8, :));
    %example = aveshape(k, :)';
   for k1 = 1 : nContrast
       sprm(k, k1, :) = fminsearchbnd(@(x) fitShiftedScaledCopy(example, x, squeeze(cprd(k, k1, :))), init, lb, ub);
       sprd(k, k1, :) = computeShiftedScaledCopy(example, squeeze(sprm(k, k1, :))); 
   end
   % compute R2
   tmp1 = sprd(k, :);
   tmp2 = cprd(k, :);
   top = sum((tmp1 - tmp2).^2);
   bottom = sum(tmp2.^2);
   r2(k) = 1-top/bottom
end


figure (2), 
for k = 1 : 3
   subplot(2, 3, k + 3), cla,  set(gca, 'colororder', copper(nContrast)), hold on
   plot(t, squeeze(sprd(k, :, :))', 'linewidth', 2), box off, ylim([0, 1])
   
end


