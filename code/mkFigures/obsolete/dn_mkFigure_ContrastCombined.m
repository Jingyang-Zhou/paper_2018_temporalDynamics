% dn_mkFigure_ContrastCombined

% DESCRIPTION -------------------------------------------------------------
% This code fit the DN model to single cells' response to a static stimulu
% but of different levels of contrast. The code also makes figure for the
% contrast part of figure 4.

%% SAVE FIGURE ?

saveFigure = 0;

%% LOAD ALBRECHT AND GEISLER DATA

dataLoc = fullfile(dn_ECoG_RootPath, 'data');
fName   = 'Figure1_ACE_Data.xlsx';
a       = xlsread(fullfile(dataLoc, fName)); % concatenated time course (over 3 cells) x contrast levels

%% VISUALIZE THE ORIGINAL DATA a

figure (100), clf, set(gca, 'colororder', copper(11)), hold on
plot(a, 'linewidth', 2), axis tight, set(gca, 'fontsize', 14)
xlabel('first dimension'), ylabel('second dimension')

%% EXTRACT AND DERIVE PARAMETERS

% MAKE STIMULUS -----------------------------------------------------------
stim = [ones(1, 200), zeros(1, 300)];
t    = 1 : 1 : length(stim);

% EXTRACT RESPONSE FOR EACH CELL ------------------------------------------
contrast_levels = a(1, 2 : end);
cellRsp{1}  = a(2 : 27, 2 : 11)./100; % the time courses have maximal repsonse 100, and we normalize to maximal response 1 here
time{1}     = a(2 : 27, 1); % in unit of millisecond                                                                                                                                           []p

cellRsp{2}  = a(30 : 50, 2 : 11)./100; % response from low to hight
time{2}     = a(30 : 50, 1);

cellRsp{3}  = a(53 : 73, 2 : 11)./100;
time{3}     = a(53 : 73, 1);

ncells      = 3;
nContrast   = length(contrast_levels);

%% INTERPOLATE THE DATA AND COMBINE THE TIME COURSES FROM THE 3 CELLS

combined_time = union(time{1}, time{2}); % the second and the third time cell have the same sampling points
combined_rsp  = [];

for icell = 1 : ncells
    for k = 1 : nContrast
        combined_rsp(icell, k, :) = interp1(time{icell}, cellRsp{icell}(:, k), combined_time, 'nearest');
    end
end

figure (100), clf
for icell = 1 : ncells
   subplot(2, 3, icell), set(gca, 'colororder', copper(nContrast)), hold on
   plot(combined_time, squeeze(combined_rsp(icell, :, :))', 'linewidth', 2), xlim([30, 150])
   title(sprintf('CELL %d', icell))
end

% Take the average across cells
mCombined_rsp = squeeze(nanmedian(combined_rsp));

subplot(2, 3, 4), set(gca, 'colororder', copper(nContrast)), hold on
plot(combined_time, mCombined_rsp, 'linewidth', 2), xlim([30, 150])

mCombined_rsp = squeeze(nanmean(combined_rsp));

subplot(2, 3, 5), set(gca, 'colororder', copper(nContrast)), hold on
plot(combined_time, mCombined_rsp, 'linewidth', 2), xlim([30, 150])
