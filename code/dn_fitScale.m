function scale = dn_fitScale(data, prm, stim, t, irfType)

% It should just be a linear regression instead of a fit

%% compute the scale for 

for k = 1 : size(data, 1)
    % compute model prediction
    [~, prd(k, :)] = dn_computeFineFit(prm(k, :), data(k, :), stim, t, irfType);
    
    % compute the scale
    scale(k) = regress(data(k, :)', prd(k, :)');
end

end