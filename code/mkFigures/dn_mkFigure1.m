% dn_mkFigure1

%% model parameteres

param = [0.1, 0, 0.1, 1, 0.2, 0, 1];
normMax_range = @(x, range) (x - mean(x(range)))./max(x - mean(x(range)));
normMax       = @(x) x./max(x);

%% summation
t    = [1 : 1000]./1000;
stim = [];

stim.summation = zeros(2, length(t)); 
stim.summation(1, 100 : 150) = 1; stim.summation(2, 100 : 200) = 1;

for k = 1 : 2
    model.summation(k, :) = dn_DNmodel(param, stim.summation(k, :), t);
end

figure (1), 
for k = 1 : 2
   subplot(3, 2, k), cla, plot(t, stim.summation(k, :), 'k-', 'linewidth', 3), hold on, 
   plot(t, model.summation(k, :), 'r-', 'linewidth', 3), box off
   set(gca, 'xtick', [0, 1], 'ylim', [0, 1], 'ytick', [0, 1])
end
%% adaptation

stim.adaptation = zeros(2, length(t));
stim.adaptation(1, [100 : 200, 300 : 400]) = 1;
stim.adaptation(2, [100 : 200, 500 : 600]) = 1;

for k = 1 : 2
    model.adaptation(k, :) = normMax_range(dn_DNmodel(param, stim.adaptation(k, :), t), 1 : 100);
end

figure (1), 
for k = 1 : 2
   subplot(3, 2, k+2), cla,  plot(t, stim.adaptation(k, :), 'k-', 'linewidth', 3), hold on, 
   plot(t, model.adaptation(k, :), 'r-', 'linewidth', 3), box off,
   set(gca, 'xtick', [0, 1], 'ylim', [0, 1], 'ytick', [0, 1])
end

%% arbitrary temporal stimuli

stim.arbitrary = zeros(2, length(t));
stim.arbitrary(1, :) = (square((t+1).*11.5)+1)/2; tmp = 1;%heaviside(t - 0.4).*0.5+0.2; 
stim.arbitrary(1, :) = stim.arbitrary(1, :).*tmp;
stim.arbitrary(1, 100 : 200) = 1;
%stim.arbitrary(2, [100 : 200, 400 : 700]) = 1; stim.arbitrary(2, [200 : 400]) = 0.5;
%stim.arbitrary(2, :) = (sin(t.*40)+1)./2;

stim.arbitrary(2, [200 : 250]) = 1;
stim.arbitrary(2, [500 : 550]) = 1;
stim.arbitrary(2, [800 : 850]) = 1;

stim.arbitrary(2, [100 : 150]) = 1;
stim.arbitrary(2, [320 : 420]) = 1;
stim.arbitrary(2, [690 : 770]) = 1;

for k = 1 : 2
    model.arbitrary(k, :) = normMax(dn_DNmodel(param, stim.arbitrary(k, :), t));
end

figure (1), 
for k = 1 : 2
   subplot(3, 2, k+4), cla,  plot(t, stim.arbitrary(k, :), 'k-', 'linewidth', 3), hold on, 
   plot(t, model.arbitrary(k, :), 'r-', 'linewidth', 3), box off,
   set(gca, 'xtick', [0, 1], 'ylim', [0, 1], 'ytick', [0, 1])
end