% sit2009_visualization

function [] = sit2009_visualize(pred, prm, sz, whichFg)

%% prep for plotting

stim = prm.stim;
rf   = prm.rf;

t_lth = size(pred{1}, 3);
mid   = round(sz/2);

max_1 = max(pred{1}(:));
max_2 = max(pred{2}(:));

pred_toplot{1} = pred{1}./max_1; % high contrast not spatially normalized
pred_toplot{3} = pred{2}./max_1; % low contrast normalized to high contrast, not spatially normalized

% for the second row
for it = 1 : t_lth
    tmp_1 = pred{1}(:, :, it); max_11 = max(tmp_1(:)); if max_11 == 0, max_11 = 1; end
    tmp_2 = pred{2}(:, :, it); max_12 = max(tmp_2(:)); if max_12 == 0, max_12 = 1; end
    
    pred_toplot{2}(:, :, it) = pred{1}(:, :, it)./max_11;
    pred_toplot{4}(:, :, it) = pred{2}(:, :, it)./max_12;
end

% for the third row


%% plotting

figure (whichFg), clf

title_txt = {'high ctr', 'high ctr norm. to sp.', 'low ctr', 'low ctr norm.to sp.+ctr.'};
% a row of movies:
for it = 1 : t_lth
    % plot stimulus
    ax1 = subplot(4, 5, 1); colormap(ax1, 'gray'), imagesc(stim.st(:, :, it), [0, 1]), axis off
    
    for k = 1 : 4
        % visualize overall response
        ax2 = subplot(4, 5, k+1); cla, colormap(ax2, 'jet'), imagesc(rf.x, rf.y, pred_toplot{k}(:, :, it), [0, 1])
        title(title_txt{k}), xlabel('(mm)'), ylabel('(mm)'), set(gca, 'fontsize', 12) 
        
        % visualize individual slice of responses (1)
        subplot(4, 5, 5 + k + 1), set(gca, 'color', 'k'), %set(gca, 'colororder', gray(sz)), 
        plot(prm.t(1 : it), squeeze(pred_toplot{k}(:, mid, 1 : it)), 'w-', 'linewidth', 2.5), hold on
        plot(prm.t(1 : it), squeeze(pred_toplot{k}(mid, mid, 1 : it)), 'r-', 'linewidth', 3), 
        plot(prm.t(1 : it), squeeze(pred_toplot{k}(end, mid, 1 : it)), 'b-', 'linewidth', 3),
        box off, xlim([prm.t(1), prm.t(t_lth)]), ylim([0, 1]), set(gca, 'ytick', ''), xlabel('time (s)')
        
        % visualize individual slice of repsonse (2)
    end
    drawnow, pause(0.1)
end

% visualize individual slice of repsonse (2)
for k = 1 : 4
   subplot(4, 5, 10 + k +1),  hold on, %set(gca, 'colororder', gray(sz/2)), hold on
   tmp = pred_toplot{k}(:, mid, :);
   for k1 = 1 : size(tmp, 1)
       tmp1 = pred_toplot{k}(k1, mid, :);
       slice_toplot(k1, :) = tmp1./max(tmp1(:));
   end
   % plot all slices
   plot(prm.t, slice_toplot(1 : mid, :), 'w-', 'linewidth', 2.5),
   % a center slice
   plot(prm.t, slice_toplot(mid, :), 'r-', 'linewidth', 3), 
   % a periphery slice
   plot(prm.t, slice_toplot(5, :), 'b-', 'linewidth', 3)
    set(gca, 'ytick', ''), xlabel('time (s)'), set(gca, 'color', 'k')
end

end