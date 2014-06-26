function [radii thetas] = polar_1D_get_radii_angles(Features, spacings)
% Get radial and angular information for the circular constructions
% associated with each feature.
%
% Arguments:
%
% Features : cell array of features
% spacings : vector of translation amounts

num_features = length(Features);
radii = zeros(num_features, 1);
thetas = zeros(num_features, 1);

for i = 1 : num_features
    F = fbshift(Features{i}, spacings(i) / 2);
    thetas(i) = 2*acos( dot(F(:,3)-F(:,2),F(:,3)-F(:,1)) / norm(F(:,3)-F(:,1)) / norm(F(:,2)-F(:,3)) );
    radii(i) = norm(F(:,2)-F(:,3)) / (2 * cos((pi-thetas(i))/2) ); % use cosine rule
end