function derivedPrm = dn_computeDerivedParams(prm, irfType)


%% compute derived parameters

t    = 0.001 : 0.001 : 10;
stim = ones(1, length(t));

% compute model responses
for k = 1 : size(prm, 1)
   [~, rsp(k, :)] = dn_computeFineFit(prm(k, :), stim, stim, t, irfType);
   
   % compute time to peak
   derivedPrm.t2pk(k) = find(rsp(k, :) == max(rsp(k, :)));
   
   % compute asymptotic response
   derivedPrm.r_asymp(k) = rsp(k, end);
end

end