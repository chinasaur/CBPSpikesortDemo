function result = CheckFeatures(Features)
% Validate a set of features e.g., if there are any NaNs then reject.
% Arguments : Features : cell array of vector/matrix features.
% Returns a boolean which is true if Features are valid.
for feature_num = 1 : length(Features)
    feature = Features{feature_num};
    if (sum(isnan(feature(:))) > 0)
        result = false;
        return;
    end
end
result = true;