function whitened_vec = WhitenVectorBlocks(vec, lens, ...
                                           num_channels, whiten_fn)

    if (size(vec, 1) ~= sum(lens) * num_channels)
        error('length of vec must be sum of lens (%d)', ...
              sum(lens) * num_channels);
    end

    vecs = vec2cell(vec, lens, num_channels);
    whitened_cells = cell(size(vecs));
    1;
    parfor example_num = 1 : length(lens)
        len = lens(example_num);
        whitened_cells{example_num} = whiten_fn(len) * vecs{example_num};    
    end
    whitened_vec = vectorize_cell(whitened_cells);
end