function [TransformParams, Magnitudes, reconstructed_data, Info] = ...
    cbp_core_wrapper(data_samples, Features, pars)

% Do many examples in parallel
%
% Arguments:
%
% data_samples : cell array of data examples
% Features : cell array of features
% pars : parameters
%
% Returns:
% TransformParams : cell array estimated transformational parameters
% Magnitudes : cell array with estimated instance Magnitudes%
% reconstructed_data : cell array reconstruction of the data examples.
% Info : information (including grid points, Dictionary, raw coeffs,
%        weights - see cbp_core2.m).

% Simply wraps around cbp_core.m
if ~isfield(pars, 'progress'), pars.progress = false; end

TransformParams = cell(length(data_samples), 1);
Magnitudes = cell(size(TransformParams));
reconstructed_data = cell(size(TransformParams));
info_flag = nargout > 3;
if (info_flag)
    %fprintf('Warning: storing info. THIS COULD BE MEMORY INTENSIVE!\n');    
    %fprintf('Press Enter to continue.\n');
    %pause;
    Info = cell(size(TransformParams));
end
true_spikes_flag = isfield(pars, 'true_spikes') && ...
                   ~isempty(pars.true_spikes);
if (isfield(pars, 'parfor_chunk_size'))
    PARFOR_CHUNK_SIZE = pars.parfor_chunk_size;
else
    PARFOR_CHUNK_SIZE = 1e3;
end
offset = 0;
num_chunks = max(1, ceil(length(data_samples) / PARFOR_CHUNK_SIZE));

% Process in chunks so as not to be so memory-intensive.
for chunk_num = 1 : num_chunks
    chunk_idx = offset + 1 : min(length(data_samples), ...
                                 offset + PARFOR_CHUNK_SIZE);
    chunk_size = length(chunk_idx);
    T = cell(chunk_size, 1);
    M = cell(chunk_size, 1);
    R = cell(chunk_size, 1);
    if info_flag
        I = cell(chunk_size, 1);
    else
        I = [];
    end
    fprintf('Processing chunk %d/%d of size %d...\n', ...
             chunk_num, num_chunks, chunk_size);

    D = data_samples(chunk_idx);
    greedy_count = false(chunk_size, 1);
    info_flag = nargout > 3;
    pfp = []; % Has to be defined for parfor...
    if pars.progress
        pfp = ParforProgress2('Running CBP...', chunk_size, 0.05);
    end
    parfor example_num = 1 : chunk_size
        tic;
        
        pars_i = pars;
        if (true_spikes_flag)
            pars_i.true_spikes = pars.true_spikes{example_num};
        end
        if (pars_i.debug_mode)
            fprintf('\n\n\tEvaluating example %d/%d\n', ...
                    example_num, chunk_size);
        end
        
        [T{example_num}, ...
         M{example_num}, ...
         R{example_num}, ...
         tmp] = ...
             cbp_core2(D{example_num}, Features, pars_i);  
         
        greedy_count(example_num) = ...
           isfield(tmp, 'greedy_soln') && ...
           tmp.greedy_soln;
        
        if (pars_i.debug_mode)
            fprintf('\n\n');
        end
        
        if (info_flag)
            I{example_num} = tmp;
        end
        
        eltime = toc;
        if pars_i.progress
            pfp.increment(example_num); 
        end
    end
    if pars.progress, delete(pfp); end
    
    % Copy them into the master data structures
    for i = 1 : chunk_size
        TransformParams{chunk_idx(i)} = T{i};
        Magnitudes{chunk_idx(i)} = M{i};
        reconstructed_data{chunk_idx(i)} = R{i};
        if info_flag  % Warning: memory intensive!
            Info{chunk_idx(i)} = I{i};
        end
    end
    fprintf('done.\n');
    fprintf('%d / %d (%0.3f percent) used greedy solution.\n', ...
            sum(greedy_count), chunk_size , 100 * mean(greedy_count));
    offset = chunk_idx(end);
end
