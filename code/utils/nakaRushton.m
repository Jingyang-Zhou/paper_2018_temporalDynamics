% naka-rushton function 

% 3 parameters:
% dm : asymptotic performance at the highest contrast
% s50: semi-saturation constant
% n : slope

function output = nakaRushton(x, t)

% x represents parameters, and t represents stimulus levels

% examples:
% t = 0 : 0.1 : 2;
% x = [];
% x.dm  = 2;
% x.s50 = 0.3;
% x.n   = 4;

output = x.dm.* (t.^x.n)./(t.^x.n + x.s50.^x.n);

%% visualize

% figure (100), clf
% semilogx(t, output), xlim([0, 1])

end