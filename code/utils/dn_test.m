% dn_test

%% USEFUL FUNCTIONS

normMax = @(x) x./max(x);

%% USE THE SAME SET OF PARAMS TO MAKE PREDICTIONS TO STIMULUS OF DIFFERENT LENGTH

t    = 0.001 : 0.001 : 1.2;
stim = zeros(1, length(t));
stim(201 : 700) = 1;

% STIMULUS OF REDUCED LENGTH:

trimRng = 1 : 700;
tTrim   = t(trimRng );
sTrim   = stim(trimRng );

%% SET UP MODEL PARAMETERS

test_prd = {};

x{1} = [0.1, 0, 0.1, 2, 0.05, 0, 1];
x{2} = [0.9, 0, 0.1, 3, 0.05, 0, 1];

for k = 1 : 2
    
    test_prd{1}(k, :) = normMax(dn_DNmodel(x{k}, stim, t));
    test_prd{2}(k, :) = normMax(dn_DNmodel(x{k}, sTrim, tTrim));
end

figure (100), clf

subplot(1, 2, 1)
plot(t, test_prd{1}(1, :)), hold on
plot(tTrim, test_prd{2}(1, :))

subplot(1, 2, 2)
plot(t, test_prd{1}(2, :)), hold on
plot(tTrim, test_prd{2}(2, :))

