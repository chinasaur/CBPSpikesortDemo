function [c_coeffs u_coeffs v_coeffs optval] = cbp_ecos_2norm(params, varargin)

opts = inputParser();
opts.addParamValue('verbose', 0);
opts.addParamValue('prec', []); % Ignored
opts.parse(varargin{:});

% Formulate problem and translate data to ECOS input format
[ecos_c, ecos_G, ecos_h, dims, ecos_A, ecos_b] = cbp2ecos_2norm(params);

% Solve!
[~, xf, ~, ~, xK] = ecos(ecos_c, ecos_G, ecos_h, dims, ecos_A, ecos_b, opts.Results);
xx = [xf; xK];
x = full(xx);
optval =  x' * [ecos_b;ecos_h];

% Translate ECOS output into problem output
ndict = size(params.dict, 2) / 3;
ldata = length(params.data);
[c_coeffs, u_coeffs, v_coeffs] = ecos2cbp_2norm(x, ndict, ldata);


function [c_coeffs, u_coeffs, v_coeffs] = ecos2cbp_2norm(x, ndict, ldata)
c_coeffs = x(1:ndict);
u0 = ndict*3 + 1 + ldata;
uind = u0 + (2:3:ndict*3);
u_coeffs = x(uind);
v_coeffs = x(uind+1);


function [ecos_c, ecos_G, ecos_h, dims, ecos_A, ecos_b] = cbp2ecos_2norm(params)
% RADII, THETA must be column vectors

dict = params.dict;
data = params.data;
ndict = size(dict, 2) / 3;
ldata = length(data);


% Cone constraints in conventional order (ecos_G in particular will have to
% be arranged to match these conventions)
dims.f = ndict;
dims.l = ndict*2;
dims.q = [ldata+1; 3.*ones(ndict,1)];


% Calculate locations for variables
slack1 = 1;
slack2 = ndict*2;
quadobj = slack2+1;
slack3 = quadobj+1;
slack4 = quadobj+ldata;
quadconstr1 = slack4+1;
quadconstr2 = slack4+ndict*3; 


ecos_c = data;
ecos_c = [ecos_c; zeros(ndict*2,1)];
ecos_h = zeros(quadconstr2,1);
ecos_h(quadobj) = 1 / sqrt(2) / params.noisesigma;


% ecos_G...

% Objective norm constraint slacks
Oncs = -speye(ldata);

% Linear constraint slacks
Lcs = -speye(ndict*2);

% Radii constraint dictionary equalities
Rcde = reshape([sparse(ldata,ndict); -dict(:,2:3:end); -dict(:,3:3:end)], ldata, ndict*3)';

% Radii constraint eye
Rceye1 = reshape([-speye(ndict);        sparse(ndict*2,ndict)],             ndict, ndict*3)';
Rceye2 = reshape([ sparse(ndict,ndict); speye(ndict); sparse(ndict,ndict)], ndict, ndict*3)';
Rceye = [Rceye1 Rceye2];

ecos_G1 = [sparse(ndict*2,ldata); sparse(1,ldata);   Oncs;                  Rcde ];
ecos_G2 = [Lcs                  ; sparse(1,ndict*2); sparse(ldata,ndict*2); Rceye];
ecos_G = [ecos_G1 ecos_G2];


% Free variables (c_coeffs)
F = sparse(-dict(:,1:3:end)');

% Radii constraints
Rc = sparse(1:ndict, 1:ndict,  params.radii);

% Linear constraints
Lc = sparse(1:ndict, 1:ndict, -params.radii.*cos(params.theta));

ecos_A = [F Rc Lc];
ecos_b = params.lambda;