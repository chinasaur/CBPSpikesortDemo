function aligned_waveform = align_peak(waveform, desired_length, mode)
% Align the peak of a waveform to the center index.
% Arguments:
% waveform : vector of waveform samples.
% desired_length : desired length of output
% mode : either {peak, centroid, biphasic} for different alignment methods.

% Returns a vector of size desired_length with peak aligned.
if(mod(size(waveform, 1), 2) == 0)
    1;
    error('waveform length must be odd!');
end

if (norm(waveform) < eps)
    aligned_waveform = waveform;    
    return;
end

if (nargin < 3)
    mode = 'peak';
end

if (desired_length > size(waveform, 1))
    waveform =  ...
        [zeros(floor((desired_length - size(waveform, 1)) / 2), ...
            size(waveform, 2)); ...
         waveform; ...
         zeros(floor((desired_length - size(waveform, 1)) / 2), ...
            size(waveform, 2))];
end

switch(mode)
    case 'peak'
        [blah max_index] = max(sum(abs(waveform), 2));
    case 'max1'
        [blah max_index] = max(max(waveform, [], 2));
    case 'min1'
        [blah max_index] = min(min(waveform, [], 2));
    case 'max'
        [blah max_index] = max(sum(waveform, 2));
    case 'min'
        [blah max_index] = min(sum(waveform, 2));
    case 'centroid'
        max_index = round(centroid(sum(abs(waveform), 2)));
    case 'biphasic'
        [blah mini] = min(sum(waveform, 2));
        [blah maxi] = max(sum(waveform, 2));
        max_index = round( (mini + maxi) / 2 );
end

aligned_waveform = zshift(waveform, ceil(size(waveform, 1) / 2) - max_index);
aligned_waveform = aligned_waveform( ...
                   (ceil(size(aligned_waveform, 1) / 2) - ...
                   floor(desired_length / 2)) : ...
                   (ceil(size(aligned_waveform, 1) / 2) + ...
                   floor(desired_length / 2)), :);

