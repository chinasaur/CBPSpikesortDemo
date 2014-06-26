function outax = sanesubplot(nrows, ncols, thisPlot, varargin)
% SANESUBPLOT   The Mathworks subplot is hard to use; this is more rational
% usage: sanesubplot(nrows, ncols, thisPlot, varargin)
% 
% Most arguments just get passed through to subplot.  The only difference is
% how thisPlot is handled.  Instead of being a counter-intuitive 1D index,
% you just use 2D indices.  If you want a subregion from rows 3:4 and cols
% 4:5, just put {3:4 4:5}.  If you want a single region, just do {x y} or 
% [x y].
%
%
% Peter H. Li 07-Nov-2011
% As required by MatLab Central FileExchange, licensed under the FreeBSD License
%

if iscell(thisPlot)
    rows = thisPlot{1};
    cols = thisPlot{2};
elseif numel(thisPlot) == 2
    rows = thisPlot(1);
    cols = thisPlot(2);
elseif numel(thisPlot) == 1
    [rows cols] = ind2sub([nrows ncols], thisPlot);
end


[R,C] = meshgrid(rows, cols);
I = sub2ind([ncols nrows], C, R);

ax = subplot(nrows, ncols, I(:), varargin{:});


if nargout > 0
    outax = ax;
end