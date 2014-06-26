function residual = ComputeWhitenedResidual(data, reconstructed_data, ...
                                            whitening_mtx_fn)

    data_vec = vectorize_cell(data);
    reconstructed_data_vec = vectorize_cell(reconstructed_data);
    residual = WhitenVectorBlocks(data_vec - reconstructed_data_vec, ...
                                  cell_length(data), size(data{1}, 2), ...
                                  whitening_mtx_fn);


end
