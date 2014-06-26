function PlotSparsityProgress(iteration_num, SparsityProgress, window_len,...
                     font_size, line_width, marker_size)

% Plot the progress of the sparsity(prior term) vs. 
% iteration in current figure.           
                    
% Arguments : 
% iteration_num : current iteration number
% SparsityProgress : vector of negative log priors vs. iteration
% window_len : sliding window length
% font_size, line_width, marker_size : plotting parameters

idx = max(1, iteration_num - window_len) : iteration_num; % iterations to plot
plot(idx,SparsityProgress(idx), '.-', ...
     'LineWidth', line_width, 'MarkerSize', marker_size);
if (iteration_num > 1)
    xlim([idx(1) idx(end)]);
end
set(gca,'FontSize', font_size);
sparsity = SparsityProgress(iteration_num);
if (iteration_num > 1)
    prev_sparsity = SparsityProgress(iteration_num - 1);
    title(sprintf('sparsity = %0.3f Per. Decr = %0.3f', sparsity, ...
          (prev_sparsity - sparsity) / abs(prev_sparsity) * 100));
else
    title(sprintf('sparsity = %0.3f',sparsity));
end
xlabel('iteration number');
ylabel('prior term');
