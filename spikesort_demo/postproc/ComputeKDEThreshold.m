function threshold = ComputeKDEThreshold(data, num_points, range, peak_width)

[bandwidth, density, xmesh, cdf] = ...
    KernelDensityEstimator(data, num_points, range(1), range(2));
%density_ddx = diff(density) ./ diff(xmesh);
% Find the last local minimum, i.e. when the first derivative is 0 and
% second derivative is positive.
%density = smooth(density, 5);
%density = density ./ sum(density) ./ mean(diff(xmesh));
if (~exist('peak_width', 'var'))
    peak_width = 21;
end
idx = (peak_width + 1) : (num_points - peak_width);
cand = true(length(idx), 1);
for i = -peak_width : peak_width
    cand = cand & density(idx) <= density(idx + i);
end
try      
    best_idx = find(cand, 1, 'last');
    threshold = xmesh(best_idx + peak_width);
catch
    fprintf('Threshold computation failed!');
    threshold = -Inf;
end

if 0
    figure;
    h = histc(data, xmesh);
    h = h ./ sum(h) ./ mean(diff(xmesh));
    plot(xmesh, [h, density], '.-', 'LineWidth', 2);
    legend('empirical', 'kde');
    if (threshold ~= -Inf)
        plot_vline(threshold, 'r-');
    end
end

if isempty(threshold), threshold = 0; end