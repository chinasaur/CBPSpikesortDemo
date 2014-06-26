function [f_value, l2_residual, sparsity] = ...
            ComputeObjectiveTerms(data, reconstructed_data, pars, ...
                                  Magnitudes)

% Compute the values of the two terms of the objective
% (the L2 data fidelity term and the negative-log-prior term)
% Arguments:
% residual_squared_norms : vector residual norms (one per example)
% pars : CBP parameters
% Magnitudes : cell array of amplitude coefficients
% prior_dbns : cell array of function handles of prior distributions for each feature's amplitudes
% Returns:
% total l2 term, total sparsity term, and then 2 vectors breaking down each
% term by example.

l2_residual = ...
    1 / (2 * pars.noise_sigma ^ 2) .* ...
        sum(ComputeWhitenedResidual(data, ...
                                reconstructed_data, ...
                                GetWhiteningMtxFn(pars)) .^ 2);
% BETA = 1e-3;
sparsity_per_example = zeros(length(Magnitudes), length(Magnitudes{1}));

% this will be num_examples x num_features
l0_lambda = -log(pars.firing_rates) .* pars.lambda;
for example_num = 1 : size(sparsity_per_example, 1)
    for feature_num = 1 : size(sparsity_per_example, 2)
        
        magnitudes = Magnitudes{example_num}{feature_num};
                
        sparsity_per_example(example_num, feature_num) = ...
            sum(magnitudes > 0) * l0_lambda(feature_num);      
    end
end
sparsity = sum(sparsity_per_example(:));
f_value = l2_residual + sparsity;

end