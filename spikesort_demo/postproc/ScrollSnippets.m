function fig = ScrollSnippets(snippets, snipcenters, varargin)

opts = inputParser();
opts.addParamValue('snipindices', 1:length(snippets));
opts.addParamValue('recons', {});
opts.KeepUnmatched = true;
opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;

data.snips   = snippets(opts.snipindices);
data.centers = snipcenters(opts.snipindices);
data.spikes  = unmatched;
if ~isempty(opts.recons)
    data.recons  = opts.recons(opts.snipindices);
else
    data.recons = [];
end

fig = gen_slider_plot('data',             data,                           ...
                      'data_length_func', @(data) length(data.snips),     ...
                      'get_data_func',    @(data, sval) data.snips{sval}, ...
                      'plot_func',        @plotsnipwrap);


                  
if nargout < 1, clear fig; end


function h = plotsnipwrap(snip, sval, data, f)
center = data.centers(sval);
opts = data.spikes;

% Precalculate how many colors to use so that we have consistent colors
% from snip to snip
opts.ncolors = 0;
if isfield(opts, 'cbp')
    opts.ncolors = max(opts.ncolors, length(opts.cbp));
end
if isfield(opts, 'clust')
    opts.ncolors = max(opts.ncolors, length(opts.clust));
end
if isfield(opts, 'true')
    opts.ncolors = max(opts.ncolors, length(opts.true));
end

% Add reconstructed snippet if available
if ~isempty(data.recons), opts.recon = data.recons{sval}; end

h = PlotSnippet(snip, center, opts);