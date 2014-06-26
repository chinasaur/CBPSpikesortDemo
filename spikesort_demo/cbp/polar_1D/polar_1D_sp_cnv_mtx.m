function CoeffMtx = polar_1D_sp_cnv_mtx(Waveforms, ...
                                        spacing, ...
                                        input_signals, ...
                                        info)

% Generate a matrix from the coefficients that, when multiplied by 
% the (vectorized) waveforms, give a reconstruction of the signal, and is 
% equivalent to multiplying the waveform's dictionary matrix with the
% vectorized coefficients.

% Arguments : 
% Waveforms : N x 1 cell array of features
% spacing : N x M spacings (M = #transformations)
% input_signals : K x 1 cell array of signals (K = # data examples)
% info : optimization info from solving coefficients. It is a K x 1 cell
%        array where each element is a struct with the fields:
%        grid_points : N x 1 cell array of gridpoints
%        raw_coeffs : raw coefficient vector for this example

% Extract the C,U,V coeffs from Xhat
num_signals = length(input_signals);
if (length(info) ~= num_signals)
    error('Info must be cell array of length = number of examples!');
end

c_coeffs = cell(num_signals, 1);
u_coeffs = cell(num_signals, 1);
v_coeffs = cell(num_signals, 1);

for signal_num = 1 : num_signals
    c_coeffs{signal_num} = info{signal_num}.raw_coeffs(1 : 3 : end);
    u_coeffs{signal_num} = info{signal_num}.raw_coeffs(2 : 3 : end);
    v_coeffs{signal_num} = info{signal_num}.raw_coeffs(3 : 3 : end);
end

waveform_total_sizes = cell_numel(Waveforms);
waveform_lens = cell_length(Waveforms);

signal_total_sizes = cell_numel(input_signals);
signal_lens = cell_length(input_signals);

shifts = round(spacing ./ 2);
num_channels = size(Waveforms{1}, 2);
num_waveforms = length(Waveforms);
num_signals = length(input_signals);

% Get the polar basis for each waveform
radii = info{1}.radii;
thetas = info{1}.thetas;

%zeros(num_waveforms, 1);
%thetas = zeros(num_waveforms, 1);
%for waveform_num = 1 : num_waveforms
%    [blah radii(waveform_num) thetas(waveform_num)] = ...
%        polar_1D_base_interp(Waveforms{waveform_num}, ...
%                             spacing(waveform_num, :));
%end

CoeffMtx = zeros(sum(signal_total_sizes), sum(waveform_total_sizes));
column_offset = 0;
1;
for waveform_num = 1 : num_waveforms
    
    % invert 3x3 polar interpolation coefficient matrix
    % We have [w- w w+] = [c u v] * CoeffMat
    CoeffMatInv = inv([1 1 1;
                       radii(waveform_num) * cos(thetas(waveform_num)), ...
                       radii(waveform_num), ...
                       radii(waveform_num)*cos(thetas(waveform_num)); ...
                       radii(waveform_num) * sin(-thetas(waveform_num)), ...
                       0, ...
                       radii(waveform_num)*sin(thetas(waveform_num))]);
    
    % Matrices that transform a waveform into (c,u,v) through 
    % left multiplication. Careful about border effects here!

    diag_idx = [-shifts(waveform_num) 0 shifts(waveform_num)];
    Cmtx = spdiags(ones(waveform_lens(waveform_num), 3) * ...
                   diag(flipud(CoeffMatInv(:, 1))), ...
                   diag_idx, ...
                   waveform_lens(waveform_num), ...
                   waveform_lens(waveform_num));    
    Umtx = spdiags(ones(waveform_lens(waveform_num), 3) * ...
                   diag(flipud(CoeffMatInv(:, 2))), ...
                   diag_idx, ...
                   waveform_lens(waveform_num), ...
                   waveform_lens(waveform_num));    
    Vmtx = spdiags(ones(waveform_lens(waveform_num), 3) * ...
                   diag(flipud(CoeffMatInv(:, 3))), ...
                   diag_idx, ...
                   waveform_lens(waveform_num), ...
                   waveform_lens(waveform_num));

    row_offset = 0;    
    for signal_num = 1 : num_signals 
        spike_idx = round(info{signal_num}.grid_points{waveform_num} + ...
                          floor(signal_lens(signal_num) / 2)) + ...
                    1; % spike indices
        if (min(spike_idx) < 1)
            fprintf('Warning: Ignoring %d spikes beyond left border\n', ...
                    sum(spike_idx < 1));
        end
        if (max(spike_idx) > signal_lens(signal_num))
            fprintf('Warning: Ignoring %d spikes beyond right border\n', ...
                    sum(spike_idx > signal_lens(signal_num)));
        end
           
        valid_idx = spike_idx > 0 & ...
                    spike_idx <= signal_lens(signal_num);        
        xlens = cell_length(info{signal_num}.grid_points);
        
        % Retrieve this waveform's (c,u,v)-coeffs for this example
        c = zeros(signal_lens(signal_num), 1);
        u = zeros(size(c));
        v = zeros(size(c));       
        
        coeff_idx = (sum(xlens(1 : waveform_num - 1)) + 1) : ...
                    sum(xlens(1 : waveform_num));
        coeff_idx = coeff_idx(valid_idx);
        spike_idx = spike_idx(valid_idx);
        1;
        c(spike_idx) = c_coeffs{signal_num}(coeff_idx);                
        u(spike_idx) = u_coeffs{signal_num}(coeff_idx);
        v(spike_idx) = v_coeffs{signal_num}(coeff_idx);
        1;
        % Create convolution matrices corresponding to these coefficients.        
        waveform_length = waveform_lens(waveform_num);        
        cnvmtx_sub_row_idx = ceil(waveform_length / 2) : ...
                             (ceil(waveform_length / 2) + ...
                             signal_lens(signal_num) - 1);                          
        C_cnvmtx = convmtx(c, waveform_length);
        CoeffCnvMtx = C_cnvmtx(cnvmtx_sub_row_idx, :) * Cmtx;
                
        U_cnvmtx = convmtx(u, waveform_length);
        CoeffCnvMtx = CoeffCnvMtx + ...
                      U_cnvmtx(cnvmtx_sub_row_idx, :) * Umtx;
        
        V_cnvmtx = convmtx(v, waveform_length);
        CoeffCnvMtx = CoeffCnvMtx + ...
                      V_cnvmtx(cnvmtx_sub_row_idx, :) * Vmtx;
        
        signal_length = signal_lens(signal_num);
        for channel_num = 1 : num_channels
            row_idx = (row_offset + ...
                       (channel_num - 1) * signal_length + ...
                       1) : ...
                      (row_offset + channel_num * signal_length);
            col_idx = (column_offset + ...
                       (channel_num - 1) * waveform_length + 1) : ...
                      (column_offset + channel_num * waveform_length);
            CoeffMtx(row_idx, col_idx) = CoeffCnvMtx;
        end
        row_offset = row_offset + signal_length * num_channels;
    end
    column_offset = column_offset + waveform_length * num_channels;    
end