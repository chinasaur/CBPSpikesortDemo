function lambda = polar_1D_lambda(firing_rates, lambda_fac)
% Choose lambda that minimizes KLdiv between
% prior distribution and Exp(lambda)
% The prior distribution for the kth cells amplitude is assumed to be:
% P(x) = (1-p) * delta(x) + p * prior_k(x)
% where p = firing_rate(k) * delta(k);

N = length(firing_rates); % number of features

if nargin < 2
    lambda_fac = 1; % factor to multiply derived lambda values by
end
if numel(lambda_fac) == 1
    lambda_fac = lambda_fac * ones(N,1);
end

lambda = zeros(N,1); % lambda values vector

for i=1:N % for each feature
    firing_prob = firing_rates(i);% * delta(i);
    %prior_amp_mean = pars.prior_amp_mean(i);
    prior_mean = firing_prob; % * prior_amp_mean;
    klfun = @(y) - (-prior_mean .* y + log(y) );
    % multiply it by user-specified factor (default=1)
    lambda(i) = lambda_fac(i) * fminbnd(klfun, 0, 1e6);    
end