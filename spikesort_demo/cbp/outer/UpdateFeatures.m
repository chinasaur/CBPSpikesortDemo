function [UpdatedFeatures step_size] = ...
    UpdateFeatures(Features, NextFeatures, step_size, pars)
% Use the step size + decay factors to update the features.
% Can optionally implement a sophisticated line search method here.

% Arguments: 
% Features : cell array of current feature estimates
% NextFeatures : cell array of optimized features.
% step_size : step size for updateing the features.
% pars : parameters (see description of outer_pars argyment in cbp_outer.m)

UpdatedFeatures = cell(size(Features));
if (CheckFeatures(NextFeatures))    
    for feature_num = 1 : length(Features)
        UpdatedFeatures{feature_num} = ...
            AdvanceWithStep(Features{feature_num}, ...
                            NextFeatures{feature_num}, ...
                            step_size);
        % Renormalize
        if (isfield(pars, 'renormalize_features') && ...
            pars.renormalize_features)
            UpdatedFeatures{feature_num} = ...
                UpdatedFeatures{feature_num} ./ ...
                norm(UpdatedFeatures{feature_num}, 'fro') .* ...
                norm(Features{feature_num}, 'fro');
                %max(abs(UpdatedFeatures{feature_num})) .* ...
                %max(abs(Features{feature_num}));
                
        end
        
        if (isfield(pars, 'step_size_decay_factor'))
            step_size = step_size * pars.step_size_decay_factor;
        end
    end
else
    fprintf('Features are invalid. Retrying.');
    UpdatedFeatures = {};    
end

end

% Advance with desired step size.
function result = AdvanceWithStep(orig, dest, step_size)
    result = orig + step_size .* (dest - orig);
end