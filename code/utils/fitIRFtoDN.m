function y = fitIRFtoDN(t,stim, params)

irf = gammaPDF_JW(t, params);

y = convCut(stim, irf, length(t));

y = y / max(y);

end
