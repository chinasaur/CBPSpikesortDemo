function PlotResidualStats(data, reconstructed_data, noise_sigma, font_size)

% Plot Residual norm distribution in current figure.

% Arguments : 
% data : cell array of signals
% reconstructed_data : cell array of reconstructed data
% noise_sigma : std. dev. of nosie
% font_size : plotting parameter

%Determine the residual norms_size (normalized by number of data points and
% noise sigma).
1;
data_vec = reshape(cell2mat(data(:)), [], 1);
reconstructed_data_vec = reshape(cell2mat(reconstructed_data(:)), [], 1);
residuals = (data_vec - reconstructed_data_vec) ./ noise_sigma;
%Rnorms = cell_fro_norms(cell_diffs(data, reconstructed_data)) .^ 2;
%residuals = Rnorms ./ cell_numel(data) ./ (noise_sigma ^ 2);    
[h c] = hist(residuals, floor(length(residuals) / 30));
h = h ./ sum(h) ./ mean(diff(c));
%bar(c, h);
plot(c, h, 'b', 'LineWidth', 2);
set(gca, 'YScale', 'log');
gauss_fit = normpdf(c, 0, 1);
hold on, plot(c, gauss_fit, 'r', 'LineWidth', 2);
if (numel(residuals) > 1)
    xlim([quantile(residuals, 1e-3), quantile(residuals, 1 - 1e-3)]);
end
try
    ylim([1e-4, max(max(h), max(gauss_fit))]);
catch
    fprintf('WARNING: error with ylim for residual dbn plot.');
end
set(gca, 'FontSize', font_size);
xlabel('residual / \sigma ');
ylabel('freq');