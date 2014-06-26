function PlotL2Progress(iteration_num, ResidualNormProgress, window_len,...
                        font_size, line_width, marker_size)
% Plot the progress of the residual norms vs. iteration in current figure.           
                    
% Arguments : 
% iteration_num : current iteration number
% ResidualNormProgress : vector of residual norms vs. iteration
% window_len : length of sliding window.
% font_size, line_width, marker_size : plotting parameters

idx = max(1, iteration_num - window_len) : iteration_num; % iterations to plot.
plot(idx, ResidualNormProgress(idx), '.-', ...
     'LineWidth', line_width, 'MarkerSize', marker_size);
if (iteration_num > 1)
    xlim([idx(1) idx(end)]);
end
set(gca, 'FontSize', font_size);
xlabel('iteration number');
ylabel('l2 term');
l2 = ResidualNormProgress(iteration_num);
if (iteration_num > 1)
    l20 = ResidualNormProgress(iteration_num - 1);
    title(sprintf('l2 = %0.3f Per. Decr = %0.3f',...
          l2, (l20 - l2) / abs(l20) * 100));
else
    title(sprintf('l2 = %0.3f', l2));
end
xlabel('L2term');
