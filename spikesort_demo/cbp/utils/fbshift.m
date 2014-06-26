% Create an m x 3 matrix where the waveform w is shifted backwards and
% forwards by n samples (center column is w itself)
function F = fbshift(w,n)

w_bwd = zshift(w, -n);
w_fwd = zshift(w, n);

F = [w_bwd(:), w(:), w_fwd(:)];

% F = [reshape([w(n+1:end,:);zeros(n,size(w,2))],numel(w),1) w(:) reshape([zeros(n,size(w,2));w(1:end-n,:)],numel(w),1)];
1;