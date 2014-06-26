function [c_coeffs u_coeffs v_coeffs cvx_optval] = cbp_cvx(params, varargin)
% Usage: [c_coeffs u_coeffs v_coeffs cvx_optval] = cbp_cvx(params, varargin)
%
% Solves CBP optimization using CVX.  Currently there are a bunch of
% commented / uncommented blocks below giving different formulations of CBP
% to feed to CVX.  In a release version, we could put these in separate
% functions and let users pick the function handle to use, or supply their
% own.

opts = inputParser();
opts.addParamValue('prec', 'default'); % From cvx_prefs
opts.KeepUnmatched = true;
opts.parse(varargin{:});
opts = opts.Results;

dict = params.dict;
data = params.data;
lambda = params.lambda;
noisesigma = params.noisesigma;
radii = params.radii;
theta = params.theta;

dim = length(lambda); % Used by CVX variable declaration

cvx_begin
  cvx_quiet(true);
  cvx_precision(opts.prec);
  
  variables c_coeffs(dim) u_coeffs(dim) v_coeffs(dim)

    % Playing with different p-norms
%   noisescale = 1 / (2 * noisesigma ^ 2);
%   modelvalue = dict * vec([c_coeffs, u_coeffs, v_coeffs]');
%   errsig = data - modelvalue;
%   errnorm = norm(errsig, 2.5);
%   errnormscale = noisescale * errnorm;
%   sparseobj = lambda' * c_coeffs;
%   minimize(errnormscale + sparseobj);
  
    % Straight 2-norm (noisescale adjusted)
%   noisescale = 1 / (sqrt(2)*noisesigma);
%   modelvalue = dict * vec([c_coeffs, u_coeffs, v_coeffs]');
%   errsig = data - modelvalue;
%   errnorm = norm(errsig);
%   sparseobj = lambda' * c_coeffs;
%   minimize(noisescale*errnorm + sparseobj);
  
    % Sum of squares - replaced with direct function
%   noisescale = 1 / (2 * noisesigma ^ 2);
%   modelvalue = dict * vec([c_coeffs, u_coeffs, v_coeffs]');
%   errsig = data - modelvalue;
%   err2norm2 = sum_square(errsig);
%   err2norm2scale = noisescale * err2norm2;
%   sparseobj = lambda' * c_coeffs;
%   minimize(err2norm2scale + sparseobj);

    % Sum of squares - original formulation broken out for profiling
%   noisescale = 1 / (2 * noisesigma ^ 2);
%   modelvalue = dict * vec([c_coeffs, u_coeffs, v_coeffs]');
%   errsig = data - modelvalue;
%   err2norm = norm(errsig, 2);
%   err2norm2 = pow_pos(err2norm, 2);
%   err2norm2scale = noisescale * err2norm2;
%   sparseobj = lambda' * c_coeffs;
%   minimize(err2norm2scale + sparseobj);
  
    % Original formulation by Chaitu
    minimize( 1 / (2 * noisesigma ^ 2) * ...
            pow_pos(norm(data - ...
            dict * ...
            vec([c_coeffs, u_coeffs, v_coeffs]'), 2), 2) + ...
            lambda' * c_coeffs )%+ const)

    % Left over simplified formulation by Chaitu
  %minimize( norm(data - ...
  %          dict * ...
  %          vec([c_coeffs, u_coeffs, v_coeffs]'), 1) + ...
  %          lambda' * c_coeffs )%+ const)
  
  subject to        

    % Cone constraint
    norms([u_coeffs v_coeffs], 2, 2) <= radii .* c_coeffs;

    % Linear constraint
    u_coeffs >= radii .* cos(theta) .* c_coeffs;

    % Amplitude constraint
%     c_coeffs <= 1.1;
cvx_end
% cvx_end runs the solver and writes c_coeffs, u_coeffs, v_coeffs, and
% cvx_optval to this workspace, from which they are returned to caller.