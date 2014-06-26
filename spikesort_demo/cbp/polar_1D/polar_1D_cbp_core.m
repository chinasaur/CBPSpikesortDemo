function [transform_params, magnitudes, info] = ...
            polar_1D_cbp_core(data, ...
                              Features, ...
                              noise_sigma, ...
                              lambda, ...
                              spacings, ...
                              info, pars)
                          

% Polar 1D CBP implementation to be used with cbp_core2.m
% Chaitanya Ekanadham 10/13/2011
%
% Arguments :
%
% data : input signal (vectorized)
% Features : cell array of num_features features
% noise_sigma : std. dev. of noise.
% lambda : weight for L1 penalty
% spacings : num_features x M matrix of transformation spacings
%            (M = #transformations)
% info (optional) : structure with any precomputed information such as:
%                   grid_points : precomputed grid points for consstructing
%                                 dictionary
%                   weights: L1 weights to be applied to each amplitude.
%                   Dictionary : precomputed dictionary
% pars : any additional parameters
%
% Returns :
% transfrom_params : estimated transformation parameters of events
% magnitiudes : estimated magnitudes "" ""
% info : optimization and dictionary information structure with fields:
%    Dictionary : interpolation dictionary of basis functions
%    grid_points : cell array of transformation param lattices with which
%                  the basis dictionary is constructed.
%    new_weights : new weights to be used in the next iteration
%                  when estimating locations/magnitudes.
%    raw_coeffs : raw coefficients from optimization that can be 
%                 left-multiplied by Dictionary to get a reconstruction of
%                 the data.

% Get the radial and angular information and store them. This will be
% used later in extracting coefficients and also in creating the
% coefficient matrix when adapting the Features themselves.
if ~exist('info', 'var') || ~isfield(info, 'radii') || ~isfield(info, 'thetas')
    if exist('pars', 'var') && isfield(pars, 'radii') && isfield(pars, 'thetas')
        info.radii = pars.radii;
        info.thetas = pars.thetas;
    else
        [info.radii, info.thetas] = polar_1D_get_radii_angles(Features, spacings);
    end
end

if exist('pars', 'var') && isfield(pars, 'precompgps') && isfield(pars, 'precompdicts')
    data_size = size(data,1);
    if length(pars.precompgps) >= data_size,
        info.grid_points = pars.precompgps{data_size};
        info.Dictionary  = pars.precompdicts{data_size};
    end
end

% Check to see if there is a precomputed dictionary. If not, compute it.
if (~exist('info', 'var') || ~isfield(info, 'Dictionary') || isempty(info.Dictionary))
    % Check to see if there are precomputed grid_points for generating the
    % dictionary. If not, compute them.
    if (~exist('info', 'var') || ...
        ~isfield(info, 'grid_points') || ...
        isempty(info.grid_points))
        info.grid_points = ...
            create_grid_points_uniform(size(data), spacings);
    end
    
    % Construct the full dictionary on the grid points.
    info.Dictionary = construct_interp_dictionary(Features, ...
        size(data), 3, @polar_1D_base_interp, spacings, @zshift, ...
        info.grid_points);        
end


% Construct vectors that specify the weights lambda, the radial info, and 
% the angular info for each triplet of coefficients.
num_blocks_per_feature = cellfun(@(C) size(C,1), info.grid_points);
lambda_vec = multirep(lambda(:), num_blocks_per_feature);
radii_vec = multirep(info.radii(:), num_blocks_per_feature);
theta_vec = multirep(info.thetas(:), num_blocks_per_feature);

% Incorporate the weights - if they're specified, they must match in dim
% to the size of lambda_vec.
if (isfield(info, 'weights') && ~isempty(info.weights))
    if (pars.debug_mode)
        fprintf('Using old weights.\n');
    end
    if (length(info.weights) ~= length(lambda_vec))
        error('Supplied number of weights (%d) must match number of blocks (%d)', ...
               length(info.weights), length(lambda_vec));
    end
    lambda_vec = info.weights;
    restrict_idx = abs(info.weights) < Inf;  
    if (sum(restrict_idx) < 2)
        if (pars.debug_mode)
            fprintf('Single coefficient: zeroing out weights!\n');
        end
        lambda_vec = zeros(size(lambda_vec));
    end
    using_old_weights = true;
else
    restrict_idx = true(length(lambda_vec), 1);    
    using_old_weights = false;
end

% If prefiltering is allowed, use simple correlation of Dictionary columns
% with data to weed out useless dictionary elements.
if (~using_old_weights && ...
    isfield(pars, 'prefilter_threshold') && ...
    pars.prefilter_threshold > 0)

    [prefilter_idx grid_point_corrs] = ...
        PrefilterDictionary(data, Features, ...
                            info.grid_points, ...
                            pars.prefilter_threshold);
    if (isfield(pars, 'debug_mode') && pars.debug_mode)
        fprintf('Prefiltering: eliminated %d of %d atom groups.\n', ...
                sum(restrict_idx & ~prefilter_idx), ...
                sum(restrict_idx));
    end    
    restrict_idx = restrict_idx & prefilter_idx;
end

if (~using_old_weights && ...
    isfield(pars, 'greedy_p_value') && pars.greedy_p_value > 0)
    % Check if we can use just one coefficient (maximally correlated one)
    % to get the residual to be very small, i.e. 
    % The squared norm of the residual should be Chi-Squared with
    % size(input, 1) degrees of freedom. If the probability of the residual
    % being greater then computed value is < greedy_p_value, then greedy is
    % not good enough. Otherwise, just optimize this coefficient!
    % For now just assume the coeff value is 1
    [best_feature_idx, best_time_idx, best_residual] = ...
        GetGreedySolutionOnLattice(Features, data, ...
                                   GetWhiteningMtxFn(pars));
    % p_value = probability that ChiSq < best_greedy_residual
    % If this is close to 1, greedy solution is bad.
    p_value = chi2cdf(best_residual ./ noise_sigma ^ 2, numel(data));
    if (pars.debug_mode)
        fprintf('Greedy solution: residual=%0.3f pvalue=%f\n', ...
                best_residual, p_value);
    end
    %fprintf('%f\n', log10(p_value));
    if (p_value < pars.greedy_p_value)
        greedy_idx = GetPolarBlockIndex(info.grid_points, ...
                                        size(data, 1), ...
                                        best_feature_idx, ...
                                        best_time_idx);                                    
        if (pars.debug_mode)
            fprintf('Restricting to 1 index: %d corresponding to', greedy_idx);
            fprintf('feature %d, time %d\n', best_feature_idx, best_time_idx);
%             figure(10), clf;
%             greedy_soln = zeros(size(data));
%             nz_idx = best_time_idx + ...
%                      (-floor(size(Features{best_feature_idx}, 1) / 2) : ...
%                        floor(size(Features{best_feature_idx}, 1) / 2));
%             sub_idx = nz_idx > 0 & nz_idx <= size(data, 1);
%             greedy_soln(nz_idx(sub_idx), :) = ...
%                 Features{best_feature_idx}(sub_idx, :);
%             plot([data, greedy_soln], '.-', 'LineWidth', 2);
%             set(gca, 'FontSize', 24);            
%             title(sprintf('p-value = %0.3f', p_value));
        end
    	restrict_idx = false(size(restrict_idx));
        restrict_idx(greedy_idx) = true;
        lambda_vec = zeros(size(lambda_vec)); % zero out all weights!
        info.greedy_soln = true;
    end
end

% Construct the "restricted" dictionary.
Dictionary_res = info.Dictionary(:, reprows(restrict_idx(:), 3));

%const = sum(-log(lambda_vec));  % add constant to make it log-likelihood
lambda_vec_res = lambda_vec(restrict_idx);
dim = length(lambda_vec_res);

WhiteningMtx = feval(GetWhiteningMtxFn(pars), numel(data));
whitened_data = WhiteningMtx * (data(:));
WhitenedDictionary_res = WhiteningMtx * Dictionary_res;

params.dict = WhitenedDictionary_res;
params.data = whitened_data;
params.lambda = lambda_vec_res;
params.noisesigma = noise_sigma;
params.radii = radii_vec(restrict_idx);
params.theta = theta_vec(restrict_idx);

%%%%%%%%%%% Solve %%%%%%%%%%%%%%%%%%
if (dim > 0)
    [c_coeffs_res u_coeffs_res v_coeffs_res info.optval] = pars.solve_fn(params, 'prec', 'low');
else
    %fprintf('All zeros. Shortcutting.\n');
    c_coeffs_res = [];
    u_coeffs_res = [];
    v_coeffs_res = [];
    info.optval = 1 / (2 * noise_sigma ^ 2) * (norm(whitened_data) ^ 2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Store raw coefficients (that can be premultiplied by dictionary to get a
% reconstruction of the data).

% If restricted idx was used, construct full coeff vectors
c_coeffs = zeros(size(restrict_idx));
c_coeffs(restrict_idx) = c_coeffs_res;
u_coeffs = zeros(size(restrict_idx));
u_coeffs(restrict_idx) = u_coeffs_res;
v_coeffs = zeros(size(restrict_idx));
v_coeffs(restrict_idx) = v_coeffs_res;

info.raw_coeffs = reshape([c_coeffs, u_coeffs, v_coeffs]', [], 1);

% Extract transformation params, magnitudes, raw coefficients.
[transform_params, magnitudes] = ...
    polar_1D_extract_fn_cvx(info.raw_coeffs, ...
                            num_blocks_per_feature, ...
                            info.Dictionary, ...
                            radii_vec, ...
                            info.grid_points, ...
                            spacings, ...
                            info.thetas, ...
                            pars.magnitude_threshold);
1;
if (pars.num_reweights > 0)
    % set the weights for the next iteration, if there is one.
    if (~isfield(pars, 'reweight_fn'))
        error('Must supply reweighting function!');
    end
    num_features = length(Features);
    if (length(pars.reweight_fn) ~= num_features)
        error('Reweight function must have %d elems.', num_features);
    end
    offset = 0;
    info.weights = zeros(size(c_coeffs));
    % EXPMIX_MEAN = 1e-3;
    for feature_num = 1 : length(Features)
        idx = offset + 1 : offset + num_blocks_per_feature(feature_num);
        info.weights(idx) = ...
            max(0, pars.reweight_fn{feature_num}(magnitudes{feature_num}));
        % Put infinite weight on coeffs already zero so they go unused in
        % next iteration
        info.weights(idx(magnitudes{feature_num} == 0)) = Inf;
        offset = idx(end);
    end
else
    offset = 0;
    info.weights = zeros(size(c_coeffs));
    for feature_num = 1 : length(Features)
        idx = offset + 1 : offset + num_blocks_per_feature(feature_num);
        % Put infinite weight on coeffs already zero so they go unused in
        % next iteration
        info.weights(idx(magnitudes{feature_num} == 0)) = Inf;
        offset = idx(end);
    end
end



%% Subroutines

function index = GetPolarBlockIndex(grid_points, data_dim, ...
    feature_index, time_index)

lens = cellfun(@(C) size(C,1), grid_points);
index = sum(lens(1 : feature_index - 1));
time_index = time_index - floor(data_dim / 2) - 1;
block_offset = find(grid_points{feature_index}(:, 1) <= time_index, ...
    1, 'last');
lhs_lattice_point = grid_points{feature_index}(block_offset);

if (block_offset < size(grid_points{feature_index}, 1) && ...
        abs(grid_points{feature_index}(block_offset + 1) - time_index) < ...
        abs(lhs_lattice_point - time_index))
    block_offset = block_offset + 1;
end

index = index + block_offset;