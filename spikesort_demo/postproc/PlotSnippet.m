function h = PlotSnippet(snippet, snipcenter, varargin)

opts = inputParser();
opts.addParamValue('dt', 1);
opts.addParamValue('startpad', 1);
opts.addParamValue('endpad', 1);
opts.addParamValue('recon', []);
opts.addParamValue('cbp',    {});
opts.addParamValue('cbpamp', {});
opts.addParamValue('cbpampthresh', []);
opts.addParamValue('clust', {});
opts.addParamValue('true',  {});
opts.addParamValue('colors', []);
opts.addParamValue('ncolors', []);
opts.addParamValue('cmf', @hsv);
opts.parse(varargin{:});
opts = opts.Results;

% Do amplitude thresholding on CBP spikes if needed
if ~isempty(opts.cbpamp) && ~isempty(opts.cbpampthresh)
    for i = 1:length(opts.cbp)
        opts.cbp{i} = opts.cbp{i}(opts.cbpamp{i} > opts.cbpampthresh(i));
    end
end


% Calculate offset
snipl = size(snippet,1);
snipmidi = (snipl-1)/2 + 1;
snipx = opts.dt*((1:snipl) - snipmidi) + snipcenter;

% Plot the data snippet
if isempty(opts.recon)
    plot(snipx, snippet);
else
    plot(snipx, snippet, 'c');
end
hold on;
axis tight

% Show reconstruction?
if ~isempty(opts.recon)
    plot(snipx, opts.recon, 'r');
end


% Note top of axis for plotting spike markers
ylim = get(gca, 'YLim');
marky = ylim(2) - 0.02*diff(ylim);

% Find CBP spikes to plot over snippet (pad a little)
findstart = snipx(1)   - opts.startpad;
findend   = snipx(end) + opts.endpad;
cbpi   = FindSpikes(findstart, findend, opts.cbp);
clusti = FindSpikes(findstart, findend, opts.clust);
truei  = FindSpikes(findstart, findend, opts.true);

% Note matching colors
if isempty(opts.colors)
    if isempty(opts.ncolors)
        opts.ncolors = max(cellfun(@length, {cbpi; clusti; truei}));
    end
    opts.colors = opts.cmf(opts.ncolors);
end

% Plot
for i = 1:length(cbpi)
    spikeindices = cbpi{i};
    if isempty(spikeindices), continue; end
    
    spiketimes = opts.cbp{i}(spikeindices);
    cbph(i) = plot(spiketimes, marky * ones(size(spiketimes)), ...
        'LineStyle', 'none', ...
        'Marker', '.', ...
        'MarkerSize', 25, ...
        'Color', opts.colors(i, :));
    if ~isempty(opts.cbpamp)
        spikeamps = opts.cbpamp{i}(spikeindices);
        for j = 1:length(spikeamps)
            text(spiketimes(j), ylim(2), sprintf('%0.2f', spikeamps(j)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'FontSize', 25);
        end
    end
end
for i = 1:length(clusti)
    spikeindices = clusti{i};
    if isempty(spikeindices), continue; end
    
    spiketimes = opts.clust{i}(spikeindices);
    clusth(i) = plot(spiketimes, marky * ones(size(spiketimes)), ...
        'LineStyle', 'none', ...
        'Marker', 'x', ...
        'MarkerSize', 20, ...
        'LineWidth', 2, ...
        'Color', opts.colors(i, :));
end
for i = 1:length(truei)
    spikeindices = truei{i};
    if isempty(spikeindices), continue; end
    
    spiketimes = opts.true{i}(spikeindices);
    trueh(i) = plot(spiketimes, marky * ones(size(spiketimes)), ...
        'LineStyle', 'none', ...
        'Marker', 'o', ...
        'MarkerSize', 15, ...
        'LineWidth', 2, ...
        'Color', opts.colors(i, :));
end

set(gca(), 'XLim', [findstart findend]);
hold off
if nargout > 0, h = gcf(); end


function ispikes = FindSpikes(snipstart, snipend, spikes)
% FindSpikes
% usage: ispikes = FindSpikes(snipcenter, snipwidth, spikes)
%
% Get indices of spikes contained within snippet at SNIPCENTER with
% SNIPWIDTH
%
ncells = length(spikes);
ispikes = cell(ncells, 1);
for i = 1:ncells
    ispikes{i} = find(spikes{i} >= snipstart & spikes{i} <= snipend);
end