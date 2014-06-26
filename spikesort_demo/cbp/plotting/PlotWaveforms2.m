function PlotWaveforms2(Features, DeltaFeaturesProgress, iteration_num, ...
                       outer_pars, font_size, line_width, marker_size)
% Plot features in INDIVIDUALLY
% The left-hand plot will plot the actual features.
% The right-hand plot will plot the % change in features vs. iteration for
% the last 10 iterations (iteration_num - 10 : iteration_num)

% Arguments : 
% Features : cell array of features
% DeltaFeaturesProgress : cell array tracking changes in Features vs.
% iteration
% iteration_num : current iteration
% outer_pars : outer parameters (see cbp_core.m)
% font_size, line_width, marker_size : plotting parameters

nrows = ceil(sqrt(length(Features)));
ncols = ceil(length(Features) / nrows);
colors = hsv(length(Features));

max_val = -Inf;
min_val = Inf;
for i = 1 : length(Features)
    max_val = max(max_val, max(Features{i}(:)));
    min_val = min(min_val, min(Features{i}(:)));
end

for i = 1 : length(Features)
    subplot(nrows, ncols, i), cla;
    xax = -floor(size(Features{i}, 1) / 2) : floor(size(Features{i}, 1) / 2);
    plot(xax, Features{i}, '.-', 'Color', colors(i, :));    
    xlim(max(cell_length(Features) ./ 2) .* [-1, 1]);
    ylim(1.25 .* [min_val max_val]);
    % If the true (golden) features are available. Plot alongside the estimates
    % in dotted lines.
    if (isfield(outer_pars, 'true_features') && ...
        ~isempty(outer_pars.true_features))
        hold on;
        plot(xax, outer_pars.true_features{i}, '.-', 'Color', colors(i, :));
    end
    set(gca, 'FontSize', font_size);
    xlabel('samples');
    if (iteration_num > 1)
        diff = DeltaFeaturesProgress(iteration_num - 1, i);
    else
        diff = nan;
    end
    title (sprintf('%d: norm=%0.3f diff=%0.3f', ...
                   i, norm(Features{i}, 'fro'), diff));
end