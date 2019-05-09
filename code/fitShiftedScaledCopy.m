% fit shifted scaled copy of the model

function lse = fitShiftedScaledCopy(example, prm, dt)

% prm
% The first entry would be scale, the second would be shift

%% make prediction for each contrast time course


input = computeShiftedScaledCopy(example, prm);


%% compare to data

lse = sum((input - dt).^2);

%% VISUALIZE

% figure (100), clf
% plot(input), hold on
% plot(dt), drawnow

end