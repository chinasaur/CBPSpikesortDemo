function VisualizeClustering(XProj, assignments, X, nchan, fig1, threshold, marker)

% XProj : snippet x PC-component matrix of projections
% assignments : vector of class assignments
% X : time x snippet-index matrix of data
if (~exist('fig1', 'var'))
    fig1 = figure;
end
if (~exist('fig1', 'var'))
    fig2 = figure;
else
  fig2 = fig1+1;
end
if (~exist('marker', 'var'))
    marker = '.';
end

N = max(assignments);
colors = hsv(N);
centroids = zeros(size(X, 1), N);
projCentroids = zeros(size(XProj,2), N);
counts = zeros(N, 1);
distances = zeros(size(X,2),1);
for i = 1 : N
  spikeIinds = find(assignments==i);
  centroids(:, i) = mean(X(:, spikeIinds), 2);
  projCentroids(:,i) = mean(XProj(spikeIinds,:)', 2);
  counts(i) = length(spikeIinds);
  distances(spikeIinds) = sqrt(sum((XProj(spikeIinds,:)'-repmat(projCentroids(:,i),1,counts(i))).^2))';
end
figure(fig1); clf;  hold on;
for i = 1 : N     %plot central cluster first
    idx = ((assignments == i) & (distances<threshold));
    plot(XProj(idx, 1), XProj(idx, 2), ...
         '.', 'Color', 0.5*colors(i, :)+0.5*[1 1 1], ...
         'Marker', marker, 'MarkerSize', 8);
end
for i = 1 : N     %plot outliers second
    idx = ((assignments == i) & (distances>=threshold));
    plot(XProj(idx, 1), XProj(idx, 2), ...
         '.', 'Color', 0.5*colors(i, :)+0.5*[1 1 1], ...
         'Marker', marker, 'MarkerSize', 8);
end
for i=1:N
  plot(projCentroids(1,i), projCentroids(2,i), 'o', 'MarkerSize', 9, 'LineWidth', 2,...
       'MarkerEdgeColor', 'black', 'MarkerFaceColor', colors(i,:));
end
xl = get(gca, 'XLim'); yl = get(gca, 'YLim');
plot([0 0], yl, '-', 'Color', 0.8 .* [1 1 1]);
plot(xl, [0 0], '-', 'Color', 0.8 .* [1 1 1]);
th=linspace(0, 2*pi, 64);
nh= plot(threshold*sin(th),threshold*cos(th), 'k', 'LineWidth', 2);
legend(nh,sprintf('spike threshold = %.1f',threshold));
axis equal
hold off
font_size = 12;
set(gca, 'FontSize', font_size);
xlabel('PC 1'); ylabel('PC 2');
title(sprintf('Clustering result (%d clusters)', N));

if (nargin < 3)
    return;
end

% Plot the time-domain snippets
figure(fig2); clf; set(gcf,'Name','Initial waveforms');
hold on;
nc = ceil(sqrt(N)); nr=ceil((N)/nc);
chOffset = 13; %**Magic vertical spacing between channels
wlen = size(X, 1) / nchan;
MAX_TO_PLOT = 1e2;
yrg = zeros(N,1);
for i = 1 : N
    Xsub = X(:, assignments == i);
    if isempty(Xsub), continue; end
    if (size(Xsub, 2) > MAX_TO_PLOT)
        Xsub = Xsub(:, randsample(size(Xsub, 2), MAX_TO_PLOT, false));
    end
    Xsub = Xsub + chOffset*floor(([1:(nchan*wlen)]'-1)/wlen)*ones(1,size(Xsub,2));
    centroid = centroids(:,i) + chOffset*floor(([1:(nchan*wlen)]'-1)/wlen);
    
    subplot(nc, nr, i);
    cla; hold on;
    plot(reshape(Xsub,wlen,[]), 'Color', 0.25*colors(i,:)+0.75*[1 1 1]);
    plot(reshape(centroid,wlen,[]), 'Color', colors(i,:), 'LineWidth', 2);
    xlim([0, wlen+1]);
    yrg(i) = max(Xsub(:))-min(Xsub(:));
    set(gca, 'FontSize', font_size);
    xlabel('time (samples)');
    title(sprintf('cell %d, snr=%.1f', i, norm(centroids(:,i))/sqrt(size(centroids,1))));
    hold off
end
% make axis ticks same size on all subplots
mxYrg = max(yrg);
for i = 1:N
  subplot(nc,nr,i);
  ylim=get(gca,'Ylim');
  ymn = mean(ylim);
  yrg = ylim(2)-ylim(1);
  set(gca,'Ylim', ymn + (ylim - ymn)*(mxYrg/yrg));
end

ip = centroids'*centroids;
dist2 = repmat(diag(ip),1,size(ip,2)) - 2*ip + repmat(diag(ip)',size(ip,1),1) +...
        diag(diag(ip));
fprintf(1,'Distances between waveforms (diagonal is norm): \n');
disp(sqrt(dist2/size(centroids,1)))
%disp(sqrt(dist2))

return

%% Subroutines

% compute expected errors for detecting signal of given
% amplitude and prior probability, assuming univariate Gaussian noise.
function [err, miss, fa] = sdtErrorRate(amp, prior)
  thresh = (log(1-prior)-log(prior) + amp^2/2)/amp;
  miss = prior * 0.5*erfc(-(thresh-amp) / sqrt(2));
  hit = prior - miss;
  fa = (1-prior) * (1 - 0.5*erfc(-thresh / sqrt(2)));
  err = 0.5*(miss/(hit + miss) + fa/(hit + fa))
  return
  
if(0)
  % simulation test:
  N=1000;
  ns = randn(round(N/prior - N), 1);
  sig = amp+randn(N, 1);
  thresh =  (log(1-prior)-log(prior) + amp^2/2)/amp;
  miss = sum(sig < thresh);
  fa = sum(ns > thresh);
  err = 0.5*(miss/N + fa/(sum(sig>=thresh)+fa))
end
  