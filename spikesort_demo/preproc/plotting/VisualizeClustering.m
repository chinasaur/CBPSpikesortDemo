function VisualizeClustering(XProj, assignments, X, nchan, fig1, fig2, marker)

% XProj : snippet x PC-component matrix of projections
% assignments : vector of class assignments
% X : time x snippet-index matrix of data
if (~exist('fig1', 'var'))
    fig1 = figure;
end
if (~exist('fig2', 'var'))
    fig2 = figure;
end
if (~exist('marker', 'var'))
    marker = '.';
end

N = max(assignments);
figure(fig1); clf;  hold on;
colors = hsv(N);
counts = zeros(N, 1);
for i = 1 : N    
    idx = assignments == i;
    counts(i) = sum(idx);
    plot(XProj(idx, 1), ...
         XProj(idx, 2), ...
         '.', 'Color', colors(i, :), ...
         'MarkerSize', 8, ...
         'Marker', marker);
     
end
font_size = 12;
set(gca, 'FontSize', font_size);
xlabel('PC 1');
ylabel('PC 2');
title(sprintf('Clustering result (%d clusters)', N));
xl = get(gca, 'XLim');
yl = get(gca, 'YLim');
plot([0 0], yl, '-', 'Color', 0.8 .* [1 1 1]);
plot(xl, [0 0], '-', 'Color', 0.8 .* [1 1 1]);
th=linspace(0, 2*pi, 64);
dof = 2;
chiMean = sqrt(2)*gamma((dof+1)/2)/gamma(dof/2);
chiVR = dof - chiMean^2;
rad= chiMean + 4*sqrt(chiVR); 
nh= plot(rad*sin(th),rad*cos(th), 'k', 'LineWidth', 2);
legend(nh,'noise, 3.9 stdevs');
axis equal
hold off

if (nargin < 3)
    return;
end

centroids = zeros(size(X, 1), N);
for i = 1 : N
    centroids(:, i) = mean(X(:, assignments == i), 2);
end

% Plot the time-domain snippets
figure(fig2); clf; hold on;
nc = ceil(sqrt(N));
%colors = hsv(nchan);
wlen = size(X, 1) / nchan;
MAX_TO_PLOT = 1e2;
for i = 1 : N
    Xsub = X(:, assignments == i);
    if isempty(Xsub), continue; end
    if (size(Xsub, 2) > MAX_TO_PLOT)
        Xsub = Xsub(:, randsample(size(Xsub, 2), MAX_TO_PLOT, false));
    end
    
    subplot(nc, nc, i);
    cla; hold on;
    offset = 0;    
    for j = 1 : nchan
        idx = offset + 1 : offset + wlen;        
        plot(Xsub(idx, :), 'Color', 0.2 .* colors(i, :) + 0.8 .* [1 1 1]);
        offset = offset + wlen;
    end
    offset = 0;
    for j = 1 : nchan
        idx = offset + 1 : offset + wlen; 
        h(j) = plot(centroids(idx, i), ...
                    'Color', colors(i, :), 'LineWidth', 2);
        lg_labels{j} = sprintf('channel %d', j);                       
        offset = offset + wlen;
    end
    xlim([0, wlen + 1]);
	ylim(1.25 .* [min(Xsub(:)), max(Xsub(:))]);
    set(gca, 'FontSize', font_size);
    xlabel('time (samples)');
    title(sprintf('Waveform %d: %d spikes', i, counts(i)));
    if ((i == 1) && (nchan > 1.5))
        l = legend(h, lg_labels);
        set(l, 'FontSize', font_size);
    end
    hold off
end

