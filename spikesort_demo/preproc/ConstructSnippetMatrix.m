function X = ConstructSnippetMatrix(data, peak_idx, pars)

% Construct a matrix where each column is a snippet of the trace centered
% about peak (vectorized across channels). Peak is computed by upsampling
% the trace, computing the root-mean-squared value across channels, and
% then smoothing. Then the trace is downsampled back to the original
% resolution.
%
% data : channel x time matrix
% peak_idx : time indices of candidate spikes.
% pars : parameters

    wlen = floor(pars.window_len / 2);
    wlen_fine = floor(pars.window_len * pars.upsample_fac / 2);
    if (pars.downsample_after_align)
        X = zeros(pars.window_len * size(data, 1), length(peak_idx));
    else
        X = zeros(pars.window_len * pars.upsample_fac * size(data, 1), ...
                  length(peak_idx));
    end
    
    % time axis
    xax = (1 : pars.window_len);
    % finely sample time axis
    fine_xax = (1 : 1 / pars.upsample_fac : pars.window_len);
    
    % Populate the matrix X one column at a time.
    parfor i = 1 : length(peak_idx)        
        % Extract snippet from the data
        x = data(:, peak_idx(i) + (-wlen : wlen));
        if pars.upsample_fac <= 1
            X(:, i) = reshape(x', [], 1);
            continue;
        end
        
        % Upsample using cubic interpolation
        x_fine = zeros(size(x, 1), length(fine_xax));
        for j = 1 : size(x, 1)
            x_fine(j, :) = interp1(xax, x(j, :), fine_xax, 'pchip');
        end
        
        % Smooth the RMS of x_fine.
        x_fine_rms = [];
        switch(pars.align_mode)
            case 'maxrms'
                x_fine_rms = sqrt(sum(x_fine .^ 2, 1));
            case 'maxabs'
                x_fine_rms = sum(abs(x_fine), 1);
            case 'max'
                x_fine_rms = sum(x_fine, 1);
            case 'min'
                x_fine_rms = sum(-x_fine, 1);
        end
        x_fine_rms = smooth(x_fine_rms, pars.smooth_len);
        
        % Align to max value of the smoothed, upsampled RMS.
        [max_val, max_idx] = max(x_fine_rms);
        % clear x_fine_rms;
        % DOWNSAMPLE BACK TO ORIGINAL RESOLUTION
        if (pars.downsample_after_align)
            new_xax = (max_idx - pars.upsample_fac * wlen) : ...
                pars.upsample_fac : ...
                (max_idx + pars.upsample_fac * wlen);
        else
            % DO NOT DOWNSAMPLE - leave in fine resolution
            new_xax = max_idx + (-wlen_fine : wlen_fine);
        end
        sub_idx = new_xax > 0 & new_xax <= size(x_fine, 2);
        lhs_pad = sum(new_xax < 1);
        rhs_pad = sum(new_xax > size(x_fine, 2));
        tmp = x_fine(:, new_xax(sub_idx));
        tmp = padarray(tmp, [0 lhs_pad], 0, 'pre');
        tmp = padarray(tmp, [0 rhs_pad], 0, 'post');
        
        % Populate matrix (vectorize across channels).
        X(:, i) = reshape(tmp', [], 1);
    end
end