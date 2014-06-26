function noise_zone_idx = GetNoiseZones(data_rms, threshold, min_len)
% Estimate noise zones from root mean square (across channels) of trace

% Identify noize zones
noise_comps = bwconncomp(data_rms < threshold);

% Determine zone lengths
comp_lens = cellfun(@length, noise_comps.PixelIdxList);

% Use only those which satisfy min length requirement.
noise_zone_idx = noise_comps.PixelIdxList(comp_lens > min_len);