function visualize_soln(TransformParams, ...
                        Magnitudes, ...
                        data, ...
                        reconstructed_data, ...
                        iteration_num, ...
                        pars)

% Function to visualize a solution to a optimization problem.
% Arguments : 
% TransformParams : num_features x 1 cell array of transformational
% parameters. For plotting purposes, the first parameter is assumed to be
% spatial (or temporal) location, and the second is assumed to be a
% dilation factor which will be reflected by the linewidth of the stem in
% the plot TODO(chaitu):  come up with better visualization!

% Color order.
colors = hsv(length(TransformParams)); % ['b' 'g' 'r' 'c' 'm' 'y ' 'k'];    

% Iteration nuber is only for setting title of plot. Optional.
if(~exist('iteration_num','var'))
    iteration_num = 0;
end

% X axis for plotting. (This is usually space or time).
x_range = floor(size(data,1) / 2) * [-1, 1];
x_axis = x_range(1) : x_range(2);

subplot(2,1,1), cla;  % Plot data and reconstruction overlayed.
plot(x_axis, data, 'b.-', 'LineWidth',2);
hold on, plot(x_axis, reshape(reconstructed_data,size(data)), ...
              'r.-', 'LineWidth', 2);

% Set fonts etc.          
set(gca,'FontSize',24);
l = legend('data','recon');
set(l,'FontSize',16);
title(sprintf('Iteration %d: err= %0.3f percent', ...
              iteration_num, ...
              mean((data(:) - reconstructed_data(:)) .^ 2) ./ ...
              var(data(:)) * 100));
xlim(x_range);

subplot(2,1,2), cla;  % Stem plot of magnitudes/transform_params
threshold = 1e-3;
MEAN_LINE_WIDTH = 2;
for feature_num = 1 : length(TransformParams)
    transform_params = TransformParams{feature_num};
    if (size(transform_params, 2) == 1)
        transform_params = [transform_params, ...
                            ones(size(transform_params, 1), 1)];
    end
    magnitudes = Magnitudes{feature_num};
    transform_params = transform_params(abs(magnitudes) > threshold, :);
    magnitudes = magnitudes(abs(magnitudes) > threshold, :);
    hold on;
    for event_num = 1 : size(transform_params, 1)
      stem(transform_params(event_num, 1), ...
           magnitudes(event_num), ...
           'Color', colors(mod(feature_num -1, length(colors)) + 1, :), ...
           'LineWidth', MEAN_LINE_WIDTH * transform_params(event_num, 2), ...
           'MarkerSize',10);
    end
    set(gca,'FontSize',24);
    xlabel('time');
    ylim([0 2]);
end

% Superimpose true spikes if available.
if (isfield(pars, 'true_spikes') && ~isempty(pars.true_spikes))
    fprintf('Superimposing true spikes.\n');
    for feature_num = 1 : length(pars.true_spikes)
        plot(pars.true_spikes{feature_num}, ...
             zeros(size(pars.true_spikes{feature_num})), ...
             'Color', colors(mod(feature_num - 1, length(colors)) + 1, :), ...
             'Marker', 'x', ...
             'LineWidth', 4, ...
             'MarkerSize', 15);
    end
end
xlim(x_range);

