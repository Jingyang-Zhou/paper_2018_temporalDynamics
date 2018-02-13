function output = dn_2Chansmodel(params, stimulus, t, dt)

% INPUTS : 
% params have fields for two beta weights, b1 and b2, and a constant
% epsilon for scaling
%% useful functions

makeIRF = @(A, B, C, t)(t/A).^8 .* exp(-t/A) - 1 / B .* (t/C).^9 .* exp(-t/C);
normSum = @(x) x./sum(x);
normMax = @(x) x./max(x);

%% defualt parameters

A = 3.29;
B = 14;
C = 3.85;

a = 2.75;
b = 11;
c = 3.18;

sigma = 0.03; % smoothing kernal
shift = 80;

%% make impulse response function

% make sustained channel impulse response
t1   = t *(1./dt);
firf = normMax(normSum(makeIRF(A, B, C, t1)));

% make transient channel impulse response
pirf = normSum(makeIRF(a, b, c, t1));
pirf = normMax(pirf - mean(pirf));

% make smoothing kernal, transofrm from neuronal to ECoG broadband response
s   = -1 : 0.001 : 1;
ker = exp(-s.^2./(2 * sigma.^2));

% compute neuronal and broadband response
for k = 1 : size(stimulus, 1)
    fcomp(k, :)  = convCut(firf, stimulus(k, :), length(stimulus));
    pcomp(k, :)  = convCut(pirf, stimulus(k, :), length(stimulus)).^2;
    % output(k, :) = normMax(params.b1 *fcomp + params.b2 * pcomp + params.epsilon);
    output(k, :) = params.b1 *fcomp + params.b2 * pcomp + params.epsilon;
    % smooth output
    output(k, :) = normMax(conv(output(k, :), ker, 'same'));
    % shift output
    output(k, :) = [zeros(1, shift), output(k, 1 : length(t) - shift)];
end


end