% dn_simulate SpikeVersusCurrent version 3

%% make Miller et al. 2009 Equation 2

% the critical "knee" frequency?
f0 = 40;      % hz
f = 1 : 1000; % frequencies

% exponents 
xl = 1.5;    % exponent in the numerator
xh = 4 - xl; % exponent in the denominator

% constant
A = 1;

% relative power spectrum:
p_f_num = A .* f.^(-xl);
p_f_den = 1 + (f./f0).^xh;

p_f = p_f_num./p_f_den;

figure (1), clf
loglog(f, p_f, 'k-', 'linewidth', 3), hold on, 
loglog([f0, f0], [10^(-9), 1], 'r--')
loglog(f, f.^(-xl), 'k:', 'linewidth', 3);
loglog(f, 9000./(f.^4), 'b:', 'linewidth', 3);
xlabel('frequency (Hz)'), ylabel('Power'), box off, ylim([10^(-9), 1])
legend('power distribution', 'knee of the power distribution', ...
    'low frequency power distribution', 'high frequency power distribution'), 
title('Component of relative power distribution')
set(gca, 'fontsize', 12)

%% make Miller et ala. 2009 figure 5

% spike rates: 
spikeArrival = poissrnd(3, [10, 100]);

figure (2), clf
subplot(2, 2, 1)
stem(spikeArrival(1, :))