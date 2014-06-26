function [precompgps precompdicts] = ...
    PrecomputeDictionaries(waveforms, spacings, snippets)

sniplens = cellfun(@(C) size(C,1), snippets);
num_chan = size(snippets{1},2);

maxsl = max(sniplens);
usls = unique(sniplens);
precompdicts = cell(maxsl,1); % Waste some pointers...
precompgps   = cell(maxsl,1);
for i = 1:length(usls)
    sniplen = usls(i);
    snipdim = [sniplen num_chan];
    
    precompgps{sniplen} = ...
        create_grid_points_uniform(snipdim, spacings);
    
    precompdicts{sniplen} = ...
        construct_interp_dictionary(waveforms, snipdim, 3, ...
        @polar_1D_base_interp, spacings, @zshift, ...
        precompgps{sniplen});
end