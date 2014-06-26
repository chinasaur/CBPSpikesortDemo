function sta = CalcSTA(data, spiketimes, window, varargin)
% Memory hungry, but vectorized for speed
% Assumes spiketimes are sorted, data is col vectors per electrode

opts = inputParser();
opts.addParamValue('dt', 1);
opts.addParamValue('dtsta', 1);
opts.addParamValue('interpmethod', 'spline');
opts.addParamValue('extrapval', 0);
opts.addParamValue('blocksize', 1000); % How many spikes to process at once; to guard against mem swap
opts.parse(varargin{:});
opts = opts.Results;

ldata = size(data,1);
nchan = size(data,2);
x = ((1:ldata) - 1) * opts.dt;

xsta = window(1):opts.dtsta:window(2);
lsta = length(xsta);
sta = zeros(lsta, nchan);

nspikes = length(spiketimes);
nblocks = ceil(nspikes / opts.blocksize);
for b = 1:nblocks
    i0 = (b-1)*opts.blocksize + 1;
    i1 = i0 + opts.blocksize;
    if i1 > nspikes, i1 = nspikes; end
    
    sts = spiketimes(i0:i1);
    xi = bsxfun(@plus, sts(:), xsta);
    snippets = interp1(x, data, xi, opts.interpmethod, opts.extrapval);
    sta = sta + reshape(sum(snippets, 1), lsta, nchan);
end

sta = sta / nspikes;