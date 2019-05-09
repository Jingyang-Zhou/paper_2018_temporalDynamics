function target = dn_fitEcog2fMRI(x, ecogrsp, mrirsp, fittype)
% function target = dn_fitEcog2fMRI(x, ecogrsp, mrirsp, fittype)
% ecogrsp:
% mrirsp:
% fittype: 'linear' or 'sqrt'

%%
switch fittype
    case 'linear'
        ecog = x.*ecogrsp;
    case 'sqrt'
        ecog = x.*sqrt(ecogrsp);
end

target = sum((ecog - mrirsp).^2);
%target = -corr(ecog', mrirsp');

%% VISUALIZE
% 
% figure (100), clf
% plot(ecog), hold on
% plot(mrirsp), drawnow

end