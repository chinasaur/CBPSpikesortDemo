function [Features, TransformParams, Magnitudes, ...
          pars outer_pars] = ...
            spikesort_outer(pars, outer_pars, data)

% Function for learning basis functions + sparse coefficients with CBP
% in the context of spike sorting. This includes re-estimation of amplitude
% prior distribution and splitting off of waveforms.

% Chaitanya Ekanadham 12/22/2011

% Arguments:

% pars : struct with fields like the following
%                   num_features : number of basis functions
%                   rates : estimated rates of the N sources (in Hz)
%                   noise_sigma : (estimated) std. dev. of the noise
%                   cbp_core_fn : function handle that solves event
%                                 transform params and magnitudes given
%                                 data and features. See cbp_core2.m for a
%                                 complete description.
%                   deltalambda_fn : function handle of the type:
%                                    [spacing lambda] = ...
%                                       deltalambda_fn(w,pars,outer_pars)
%                                    that computes spacings and lambda 
%                                    values for a set of waveforms and
%                                    parameters
% outer_pars : struct with fields:
%                num_iterations : number of outer iterations
%                batch_size : number of examples to optimize per iteration
%                             (default = 1)
%                init_features : initial estimate of basis functions
%                CoeffMtx_fn : function handle of the type:
%                          S = CoeffMtx_fn(Features, spacing, data, ...
%                                          optimization_info)
%                          that generates a matrix S that can be multiplied
%                          by the vectorized features to get a
%                          reconstruction of the data using the solved CBP
%                          coefficients. This is used to optimize the
%                          features holding the coefficients fixed.
%                adjust_wfs : (OPTIONAL) function that resizes/aligns 
%                             waveforms appropriately
%
% data : cell array of data. At each iteration, outer_pars.batch_size
%        examples, each an n x 1 vector, will be chosen randomly from this
%        set
%
% Returns:
% Features : cell array of optimized features (waveforms)
% TransformParams : cell array estimated transformational parameters
% Magnitudes : cell array with estimated instance Magnitudes

% Initialize the basis to specified initial point
Features = outer_pars.init_features;
num_features = length(Features);

% Things to track
Featuresprogress = cell(outer_pars.num_iterations, 1);
DeltaFeaturesProgress = zeros(outer_pars.num_iterations, pars.num_features);
ResidualNormProgress = zeros(outer_pars.num_iterations, 1); % L-2 term
SparsityProgress = zeros(outer_pars.num_iterations, 1); % sparsity term
Fprogress  = zeros(outer_pars.num_iterations, 1); % overall objective

if (isfield(outer_pars, 'step_size'))
    step_size = outer_pars.step_size;
else
    step_size = 1;
end

% Workaround incredibly annoying slow figures (at least in Linux)
if ~ishandle(1), figure(1); end
set(1, 'Position', [974 200 943 368]);
if ~ishandle(2), figure(2); end
set(2, 'Position', [974 798 273 302]);
if ~ishandle(3), figure(3); end
set(3, 'Position', [1255 798 662 300]);
if ~ishandle(4), figure(4); end
set(4, 'Position', [543 350 423 750]);
if ~ishandle(5), figure(5); end
set(5, 'Position', [5 350 530 748]);

1;

for iteration_num  = 1 : outer_pars.num_iterations
    pars.spacings = polar_1D_delta(Features, pars.accuracy);
    
    % Update progress variables to track performance.    
    if 1%(isfield(pars, 'debug_mode') && pars.debug_mode)
        fprintf('Spacings: %s Lambdas: %s\n',mat2str(pars.spacings, 3), ...
                mat2str(pars.lambda, 3));    
    end
    if (min(pars.spacings(:)) <= 0 || min(pars.lambda) <= 0)
        error('Invalid delta/lambda values!');
    end
    
    % Randomly sample the data
    data_sample = data(randsample(length(data), outer_pars.batch_size)); 

    % Solve CBP
    tic;
    [TransformParams, Magnitudes, reconstructed_data, info] = ...
        cbp_core_wrapper(data_sample, Features, pars);
    coeff_time = toc;
    Featuresprogress{iteration_num} = Features;  % current feature
    [Fprogress(iteration_num), ...
     ResidualNormProgress(iteration_num), ...
     SparsityProgress(iteration_num)] = ...     
        ComputeObjectiveTerms(data_sample, reconstructed_data, ...
                              pars, Magnitudes);    

    % Print some iteration information.
    fprintf('Iteration %d: Solved coeffs:\tf = %0.3f\tl2 = %0.3f\ts = %0.3f\tTime=%0.3f\n', ...
            iteration_num, ...
            Fprogress(iteration_num), ...
            ResidualNormProgress(iteration_num), ...
            SparsityProgress(iteration_num), ...
            coeff_time);        
           
    % If objective increased, return out.
    if (isfield(outer_pars, 'stop_on_increase') && ...
        outer_pars.stop_on_increase && ...
        iteration_num > 1 && ...
        Fprogress(iteration_num) > Fprogress(iteration_num - 1))
        fprintf('Objective value increased unexpectedly. Exiting.');
        return;
    end       

    % Plot the status
    if (isfield(outer_pars, 'plotevery') && ...
        (mod(iteration_num -  1, outer_pars.plotevery) == 0))  % Plot        
        plot_cbp_status(iteration_num, Features, DeltaFeaturesProgress, ...
                        ResidualNormProgress, SparsityProgress, data_sample, ...
                        reconstructed_data, TransformParams, Magnitudes, ...
                        outer_pars, pars);
    end            
        
    % Create the coefficient matrix   
    tic;
    CoeffMtx = outer_pars.CoeffMtx_fn(Features, pars.spacings, data_sample, ...
                                      info);                                      
    mtx_time = toc;
    fprintf('Generated coeff mtx in %0.3f seconds\n', mtx_time);
    if (isfield(outer_pars, 'check_coeff_mtx') && ...
        outer_pars.check_coeff_mtx)
       % Check for consistency of CoeffMtx.
       reconstructed_data_vec = vectorize_cell(reconstructed_data);
       reconstructed_data_vec2 = CoeffMtx * vectorize_cell(Features);
       if (norm(reconstructed_data_vec - reconstructed_data_vec2, 2) > ...
               1e-6)
           1;
            %error('There is an inconsistency with coeff matrix.');
       end                      
    end                  
    
    % Zero out any small values in the matrix - this will greatly increase
    % speed for basis solving.
    if (isfield(outer_pars, 'coeff_mtx_threshold') && ...
        outer_pars.coeff_mtx_threshold > 0)
        small_idx = abs(CoeffMtx) < outer_pars.coeff_mtx_threshold;        
        CoeffMtx(small_idx) = 0;
        fprintf('Zeroed out %d values\n', sum(small_idx(:)));
    end
    
    % Solve the basis functions        
    tic;
    %[success, NextFeatures, reconstructed_data] = ...
    %    solve_basis(data_sample, CoeffMtx, cell_fro_norms(Features), ...
    %                Features, pars.noise_sigma);
    [success, NextFeatures, reconstructed_data] = ...
        solve_basis_grad_step(data_sample, CoeffMtx, ...
                              Features, pars, outer_pars);                                
    basis_time = toc;
    if (~success)
        error('There was an error in the basis optimization. This could be because some features coeffs are 0 due to large lambda.');
    end
    
    % Update the change-in-waveform vector
    DeltaFeaturesProgress(iteration_num, :) = cell_mse(Features, NextFeatures)';
    
    % Update the features
    OldFeatures = Features;
    [Features step_size] = UpdateFeatures(Features, NextFeatures, ...
                                          step_size, outer_pars);
                          
    if (isempty(Features))
        fprintf('Features are invalid. Retrying.');        
        Features = OldFeatures;
        continue;
    end
    
    % If step size ~= 1, we need to recalculate the reconstructed data.    
    reconstructed_data = vec2cell(CoeffMtx * vectorize_cell(Features), ...
                                  cell_length(data_sample), ...
                                  size(data_sample{1}, 2));    
                                                                                          
    [new_obj_value, ...
     ResidualNormProgress(iteration_num), ...
     SparsityProgress(iteration_num)] = ...
        ComputeObjectiveTerms(data_sample, reconstructed_data, ...
                              pars, Magnitudes);

    fprintf('Iteration %d: Solved basis:\tf = %0.3f\tl2 = %0.3f\ts = %0.3f\tTime=%0.3f\n', ...
            iteration_num, ...
            new_obj_value, ...
            ResidualNormProgress(iteration_num), ...
            SparsityProgress(iteration_num), ...
            basis_time);

    % Terminate if objective value has increased.
    if (isfield(outer_pars, 'stop_on_increase') && ...
        outer_pars.stop_on_increase && ...
        new_obj_value > Fprogress(iteration_num))
        fprintf('Objective value increased unexpectedly. Exiting.');
        return;
    end        
    Fprogress(iteration_num) = new_obj_value;         

    % Adjust the lengths of the features    
    if(isfield(outer_pars, 'adjust_wfsize'))
        Features = outer_pars.adjust_wfsize(Features);             
    end    
end

end