function PlotMagnitudeStats(Magnitudes)

% Plot Magnitude stats in current figure
% Arguments : 
% Magnitudes : cell array where each element is a num_features x 1
% cell array of magnitudes
% PriorParams : structure containing the mixture prior paramters.
% font_size, line_width : plotting parameters

% Collect amplitudes by feature
MagnitudesPerFeature = hist_coeffs(Magnitudes);
num_features = length(MagnitudesPerFeature);
%num_cols = min(num_features, 3);
%num_rows = ceil(num_features / num_cols);
num_rows = ceil(sqrt(num_features));
num_cols = ceil(num_features / num_rows);

for feature_num = 1 : length(MagnitudesPerFeature)
    subplot(num_rows, num_cols, feature_num), cla;
    % Histogram amplitudes
    mags = MagnitudesPerFeature{feature_num};
    mags = mags(mags > 0);
    if (isempty(mags))
        continue;
    end
    [h, c] = hist(mags, 50);%floor(length(mags) / 30));
    if (isempty(c) || isempty(h))
        continue;
    end
    h = h ./ sum(h) ./ mean(diff(c));
    bar(c, h); 
end