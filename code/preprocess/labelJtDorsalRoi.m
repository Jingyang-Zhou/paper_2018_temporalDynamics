% label jt electrodes

function [index, dorsalLabels] = labelJtDorsalRoi()

names = {'V3AB', 'IPS', 'anterior IPS'};

nElectrodes    = 118;
labels         = cell(1, nElectrodes); 
ipsElectrodes  = [46, 47, 50, 55, 54];
v3abElectrodes = [49, 53, 91, 116, 117];


%% create label

for k = 1 : nElectrodes
   if (k >= 33) & (k <= 45)
       labels{k} = names{3};
   elseif sum(ismember(k, ipsElectrodes)),
       labels{k} = names{2};
   elseif sum(ismember(k, v3abElectrodes)),
       labels{k} = names{1};
   end
end

%% find the electrodes that were labeled

index = find(~cellfun(@isempty, labels));
dorsalLabels = labels(index);

end