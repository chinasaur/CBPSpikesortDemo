function norm = windowed_norm(data, windowlength, step)
% WINDOWED_NORM
% usage: norm = windowed_norm(data, windowlength, [step])
%
% Assumes data is NCHANNELSxDATALEN
% Default step size is WINDOWLENGTH / 8
%

if size(data,2) < windowlength,
    norm = zeros(size(data,1), 0);
    return
end

if nargin < 3
    step = round(windowlength / 8);
end


% Square and cumulative sum along row(s)
s = cumsum(data.^2, 2); 

% Take first sample on the front that has full WINDOWLENGTH
ss = s(:,windowlength);

% Subtract to get sumsquare within sliding window, spaced according to
% STEP, append to result
ss1 = s(:, step+windowlength:step:end);
ss2 = s(:, step             :step:end-windowlength);
ss = [ss (ss1 - ss2)];

norm = sqrt(ss);