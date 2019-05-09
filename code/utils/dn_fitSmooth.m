function target = dn_fitSmooth(sigma, data, model)
% function target = dn_fitSmooth(sigma, data, model)

normalize = @(x) (x - mean(x([1 : 200, 1000 : end])))./max(x - mean(x([1 : 200, 1000 : end])));
s      = -1 : 0.001 : 1;
ker    = exp(-s.^2./(2 * sigma.^2));
nirf   = normalize(conv(model, ker, 'same'));
target = sum((nirf - data).^2);

% figure (100), clf
% plot(nirf), hold on
% plot(data), 
% plot(ker), drawnow
end