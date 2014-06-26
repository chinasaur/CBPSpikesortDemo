function [ecos_c, ecos_G, ecos_h, dims, ecos_A, ecos_b] = cbp2ecos_2norm_cconst(dict, data, lambda, noisesigma, radii, theta, cconst)
% RADII, THETA must be column vectors
% I don't think this is finished; use QCML instead

ndict = size(dict, 2) / 3;
ldata = length(data);


% Cone constraints in conventional order (ecos_G in particular will have to
% be arranged to match these conventions)
dims.f = 0;
dims.l = ndict*3;
dims.q = [ldata+1; 3.*ones(ndict,1)];

ecos_c = data - sum(dict(:,1:3:end)).*cconst;
ecos_c = [ecos_c; radii.*cconst; -radii.*cos(theta).*cconst];

%% Rest still TODO

% Calculate locations for variables
slack1 = 1;
slack2 = ndict*2;
quadobj = slack2+1;
slack3 = quadobj+1;
slack4 = quadobj+ldata;
quadconstr1 = slack4+1;
quadconstr2 = slack4+ndict*3; 



ecos_h = zeros(quadconstr2,1);
ecos_h(quadobj) = 1 / sqrt(2) / noisesigma;


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
Rc = sparse(1:ndict, 1:ndict,  radii);

% Linear constraints
Lc = sparse(1:ndict, 1:ndict, -radii.*cos(theta));

ecos_A = [F Rc Lc];
ecos_b = lambda;
