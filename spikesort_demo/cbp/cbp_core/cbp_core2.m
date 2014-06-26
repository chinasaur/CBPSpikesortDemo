function [transform_params, ...
          magnitudes, ...
          reconstructed_signal, ...
          optimization_info] = cbp_core2(input_signal, Features, pars)

% Core CBP function: Version 2.0
% Chaitanya Ekanadham 10/14/2011

% Arguments

% input_signal : input signal (vectorized)
% Features : cell array of base features
% pars : parameter structure with fields:
%           lambda : lambda values
%           spacings : N x M spacings (M=# transforms) used in interpolation
%           cbp_core_fn : user-specified function of the form:
%            [transform_params, magnitudes, raw_coeffs, info] = ...
%                     cbp_core_fn(data, ...
%                                 Features, ...
%                                 noise_sigma, ...
%                                 lambda, ...
%                                 spacings, ...
%                                 info)
%             info should be a struct containing the Dictionary, raw
%             coeffs, grid points and weights.
%
%             (See polar_1D_cbp_core.m for an example).
%
%           num_reweights : number of reweights
%           reweight_fn : reweighting function for IRL1 (see line 23)
%
% Returns:
%
% TransformParams : estimated transformational parameters
% Magnitudes : estimated instance Magnitudes
% reconstructed_data : reconstruction of the data example.
% optimization_info : info about the optimization.
%                     Should have fields:
%        grid_points : cell array of grid points in transformation space
%                      at which the dictionary was sampled.
%
%        Dictionary : dictionary of basis fns used in this optimization
%
%        raw_coeffs : raw coefficients optimized for this example. When
%        left-multiplied by the Dictionary, this should yield a
%        reconstruction of the data.
%
%        new_weights : weights to be used in the next iteration.

counter = 0;

% Initialize the optimization_info and lambda
optimization_info = [];
lambda = pars.lambda;

if (~isfield(pars, 'num_reweights'))
    pars.num_reweights = 0;
end
magnitudes = {};
if (~isfield(pars, 'per_change_tol'))
    pars.per_change_tol = 1e-6;
end

while (counter <= pars.num_reweights)  % Main iterative reweighting loop
    
    % Estimate transform params/magnitudes.
    old_mags = cell2mat(magnitudes(:));
    [transform_params, magnitudes, optimization_info] = ...
        pars.cbp_core_fn(input_signal, Features, pars.noise_sigma, ...
                         lambda, pars.spacings, optimization_info, pars);
    
    % Reconstruct the data using the dictionary and coefficients.
    reconstructed_signal = reshape(optimization_info.Dictionary * ...
                                   optimization_info.raw_coeffs, ...
                                   size(input_signal));
    
    % Calculate the percent change in the solution for this iteration
    mvec = cell2mat(magnitudes(:));
    if (~isempty(old_mags))          
       per_change = mean((old_mags - mvec).^2) ./ var(old_mags) * 100;
    else
       per_change = 100;
    end
    
    % Number of nonzero coefficients in this iteration
	num_coeffs = sum(mvec > 0);
    
    % Debugging : Plot the current solution
    if (isfield(pars, 'debug_mode') && pars.debug_mode)
      visualize_soln(transform_params, magnitudes, ...
                     input_signal, reconstructed_signal, counter, pars);
      fprintf('Log10 percent change: %0.3f, # coeffs: %d\n', ...
              log10(per_change), num_coeffs);
      pause(0.25);
    end
    
    % If solution is not changing or there are no nonzero coeffs, break.
    if (per_change < pars.per_change_tol || num_coeffs == 0)
        break;
    end
    counter = counter + 1;
end

% Re-optimize remaining coefficients with no sparsity penalty!
if (isfield(pars, 'debug_mode') && pars.debug_mode)
    fprintf('Reoptimizing %d/%d nonzero coefficients.\n', ...
            sum(optimization_info.weights(:) < Inf), ...
            numel(optimization_info.weights));
end
optimization_info.weights(optimization_info.weights < Inf) = 0;
1;
[transform_params, magnitudes, optimization_info] = ...
    pars.cbp_core_fn(input_signal, Features, pars.noise_sigma, ...
                     lambda, pars.spacings, optimization_info, pars);
% Reconstruct the data using the dictionary and coefficients.
reconstructed_signal = reshape(optimization_info.Dictionary * ...
                               optimization_info.raw_coeffs, ...
                               size(input_signal));    

% Do a greedy solution and compare this with the optimal value.                           
if (isfield(pars, 'compare_greedy') && pars.compare_greedy)
    [transform_params, magnitudes, ...
	 reconstructed_signal, optimization_info] = ...
        CompareWithGreedySolution(input_signal, ...
                                  Features, ...
                                  lambda, ...
                                  transform_params, ...
                                  magnitudes, ...
                                  reconstructed_signal, ...
                                  optimization_info, ...
                                  pars);
end
                                                     
% Debugging : Plot final solution
if (isfield(pars, 'debug_mode') && pars.debug_mode)
    fprintf('Plotting final solution.\n');
    visualize_soln(transform_params, magnitudes, ...
                   input_signal, reconstructed_signal, Inf, pars);                
    pause;
end
