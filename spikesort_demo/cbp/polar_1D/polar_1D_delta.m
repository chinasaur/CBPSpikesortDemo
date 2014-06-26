function [delta lambda_vals] = polar_1D_delta(W, epsilon)
% Compute the best delta to use for a given set of waveforms, given a 
% desired accuracy

MIN_DELTA = 1; % minimum half-spacing
MAX_DELTA = 25; % maximum half-spacing

N = length(W); % number of features
delta = zeros(N,1); % spacings vector

for i=1:N % for each feature
    dx = MIN_DELTA;    
    while (dx <= MAX_DELTA) % increase spacing until error exceeds desired accuracy
        
        [CUV r theta] = polar_1D_base_interp(W{i},2*dx); % CUV will be 3D. Each slice same size as W{i}
        
        % These are the true shifts (vectorized if multi-dim)
        Atrue = zeros(numel(W{i}),dx+1);
        for j=1:size(Atrue,2)
            Atrue(:,j) =reshape(zshift(W{i},j-1),numel(W{i}),1); % shift by amounts 0,1,2,...,dx
        end
        
        % These are the interpolated shifts (vectorized if multi-dim)
        thvec = linspace(0,theta,size(Atrue,2));
        Ahat = reshape(CUV,numel(W{i}),3)*[ones(1,size(Atrue,2));...
                    r*cos(thvec);...
                    r*sin(thvec)];
        
        % Compute error (maximum percent difference taken across intermediate integer shifts)
        err = max(mean((Atrue-Ahat).^2)./var(Atrue,[],1)*100);
        if ( err > epsilon)
            dx = dx-1; % back up and break out if error exceeds desired accuracy
            break;
        end
        dx = dx + 1;                
    end
    delta(i) = 2 * dx;
end