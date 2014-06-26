function [best_feature_idx, best_time_idx, best_residual] = ...
    GetGreedySolutionOnLattice(Features, data, whiten_mtx_fn)

% Shift Features on the lattice and figure out the best possible residual
% by taking the data and subtracting a shifted version of each feature.
% Return the best feature index, time index, and residual value (squared L2
% norm of the residual).

num_features = length(Features);
best_feature_idx = -1;
best_time_idx = -1;
best_residual = Inf;
%WhitenMtx = whiten_mtx_fn(size(data, 1));
for feature_num = 1 : num_features
    ShiftMtx = zeros(size(data, 1), size(data, 2), size(data, 1));

    for channel_num = 1 : size(data, 2)        
        A = convmtx(Features{feature_num}(:, channel_num), size(data, 1));
        
        ShiftMtx(:, channel_num, :) = ...
            reshape( ...
              A(ceil(size(Features{feature_num}, 1) / 2) : ...
              ceil(size(Features{feature_num}, 1) / 2) + ...
              size(data, 1) - 1, :), size(data, 1), 1, size(data, 1));
    end
    residuals = (repmat(data, [1, 1, size(ShiftMtx, 3)]) - ShiftMtx);
    whitened_residuals = residuals;
    %whitened_residuals = WhitenMtx * residuals;    
    criteria = sum(squeeze(sum(whitened_residuals .^ 2, 2)), 1)';
    if (min(criteria) < best_residual)
        [best_residual, best_time_idx] = min(criteria);
        best_feature_idx = feature_num;
    end
end
    
    
    
    