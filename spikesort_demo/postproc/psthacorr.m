function [dts inds] = psthacorr(spiketimes, t0, t1, varargin)
if nargin < 2, t0 = -0.001; t1 = 0.03; end
if nargin < 4, show = nargout < 1; end
dropzero = true;
psthxcorr(spiketimes, spiketimes, t0, t1, dropzero, varargin{:});