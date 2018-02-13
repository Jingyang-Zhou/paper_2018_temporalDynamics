function Vx = sit2009_normalization(Ix, t, stim_t, g0, C, rf_lin_sp, rf_norm_sp, stim_case)
%% test:


%% model parameters

% I:  stimulus
% V:  the output response
% A:  linear repsonse, a convolution between the stimulus and the linear pool (numerator)
% B:  linear response, a convolution between the stimulus and the
%     normalization pool
% g0: a constant
% C:  a constant

% stim_case = 'offset'; % 'onset', 'steadystate', 'offset'

%% model computation

Ax = []; Bx = []; Vx = []; t_lth = length(t);
% compute spatial (static) output
for it = 1 : t_lth
    Ax(:, :, it) = conv2(rf_lin_sp, squeeze(Ix(:, :, it)), 'same');  % A(x)
    Bx(:, :, it) = conv2(rf_norm_sp, squeeze(Ix(:, :, it)), 'same'); % B(x)
end

switch stim_case
    case 'onset' % response when stimulus goes from 0 to 1
        for it = 1 : t_lth
            component_1  = Ax(:, :, it)./(g0*(1 + Bx(:, :, it)));
            component_2  = -g0/C .* (1 + Bx(:, :, it));
            Vx(:, :, it) = component_1.*(1 - exp(component_2 * t(it)));
        end
        
    case 'steadystate' % static state response to a lasting stimulus
        for it = 1 : t_lth
            component1 = Ax(:, :, it);
            component2 = g0.*(1 + Bx(:, :, it));
            Vx(:, :, it) = component1./component2;
        end 
    case 'offset' % response when stimulus goes from 1 to 0
        t0_idx = find(stim_t == 0, 1)-1;
        t0     = t(t0_idx);
        
        % compute (steady state response) V0:
        denominator = Ax(:, :, t0_idx).*Ix(:, :, t0_idx);
        numerator   = g0.*(1 + Bx(:, :, t0_idx));
        Vx(:, :, 1 : t0_idx) = repmat(numerator./denominator, [1, 1, t0_idx]);
        
        for it = t0_idx + 1 : t_lth
            component1 = 90/C*(t(it) - t0);
            Vx(:, :, it) = Vx(:, :, t0_idx).*exp(-component1);
        end
end

%% Visualize

figure (1), clf
subplot(1, 2, 1)
for it = 1 : t_lth
    imagesc(Bx(:, :, it), [0, 1]), pause(0.05)
end

subplot(1, 2, 2)
plot(t, squeeze(Vx(15, :, :))'), xlabel('time (s)')
