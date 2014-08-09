function [spike_times, spike_amps, recon_snippets] = ...
    SpikesortCBP(snippets, ...
                 snippet_centers, ...
		 waveforms, ...
                 cbp_outer_pars, ...
                 cbp_pars)
%
%function [waveforms, spike_times, spike_amps, recon_snippets] = ...
%    SpikesortCBP(snippets, ...
%                 snippet_centers, ...
%                 learn_flag, ...
%                 cbp_outer_pars, ...
%                 cbp_pars)
%
% Do spike sorting with CBP. 
% snippets : cell array of time x chan voltage trace snippets
% snippet_centers : vector of time indices into original trace associated
%                   with the center of each snippet 
% learn_flag : whether or not to learn waveform shapes
% cbp_outer_pars : pars for learning waveforms
% cbp_pars : pars for CBP inference.
%
% Returns: 
% spike_times : cell array of spike times
% spike_amps : cell array of spike magnitudes
% recon_snippets : cell array of reconstructed snippets

starttime = tic;

% Precompute dictionary, radii, thetas
cbp_pars.spacings = ...
    polar_1D_delta(waveforms, cbp_pars.accuracy);

[cbp_pars.radii, cbp_pars.thetas] = ...
    polar_1D_get_radii_angles(waveforms, cbp_pars.spacings);

[cbp_pars.precompgps cbp_pars.precompdicts] = ...
    PrecomputeDictionaries(waveforms, cbp_pars.spacings, snippets);

if (numel(cbp_pars.lambda)==1)
  cbp_pars.lambda = cbp_pars.lambda * ones(numel(waveforms), 1);
end

if (numel(cbp_pars.reweight_fn)==1)
  fn = cbp_pars.reweight_fn;
  cbp_pars.reweight_fn = cell(numel(waveforms),1);  
  for n=1:numel(waveforms)
    cbp_pars.reweight_fn{n} = fn;
  end
end

fprintf('Inferring spikes for whole data set...\n');
fprintf('CBP parameters: Spacing=%s Lambda=%s\n', ...
    mat2str(cbp_pars.spacings, 3), mat2str(cbp_pars.lambda, 3));
[spike_times_cell, spike_amps_cell, recon_snippets] = ...
    cbp_core_wrapper(snippets, waveforms, cbp_pars);
clear all_info; % save memory

% Convert cell arrays into vectors
[spike_times, spike_amps] = ConvertSpikeTimesFromCell(spike_times_cell, ... 
    spike_amps_cell, snippet_centers);

toc(starttime);

