function [wf perm] = OrderWaveformsByNorm(wf)
% ORDERWAVEFORMSBYNORM
% usage: [wf perm] = OrderWaveformsByNorm(wf)
%
% Input WF should be an MxN matrix of input waveforms, where M is the
% length of the waveforms and N is the number of waveforms.
%
% Returns the permuted waveforms matrix and the permutation used to put
% them in 2norm order.
%

nwf = size(wf,2);
norms = zeros(1,nwf);
for i = 1:nwf
    norms(i) = norm(wf(:,i));
end

[dummy, perm] = sort(norms);
wf = wf(:,perm);