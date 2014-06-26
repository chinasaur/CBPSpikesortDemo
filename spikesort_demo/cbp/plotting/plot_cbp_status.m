function plot_cbp_status(iteration_num, Features, DeltaFeaturesProgress, ...
                         ResidualNormProgress, SparsityProgress, ...
                         data, reconstructed_data, ...
                         TransformParams, Magnitudes, outer_pars, pars)
                          
% Plot the optimization status
% Chaitanya Ekanadham 06/07/2011
%
% Arguments:
% iteration_num : iteration number
% Features : current waveforms
% DeltaFeaturesProgress : percent changes in waveforms as a fn of iteration
% ResidualNormProgress : reconstruction error as a function of iteration 
% SparsityProgress : value of the prior term as a function of iteration
% data : current data_sample (cell array)
% reconstructed_data : CBP's reconstruction of the data (same size a data).
% TransformParams : current transformation parameters of estimated events
% Magnitudes : current magnitudes of estimated events in data
% outer_pars : outer parameters (see cbp_core.m)
% pars : inner parameters (see cbp_core.m)
% info : optimization info (optional)

font_size = 12; % fontsize
line_width  = 2; % line width
marker_size = 20; % marker size

% Figure 1 : Plot waveforms and percent change vs. iteration and F progress
sfigure(1); clf;
if 1 %(length(Features) < 5)
    PlotWaveforms2(Features, DeltaFeaturesProgress, iteration_num, ...
                  outer_pars, font_size, line_width, marker_size);
else
    if (isfield(outer_pars, 'true_features'))
        SubplotMultiChannelFeatures(Features, outer_pars.true_features);
    else
        SubplotMultiChannelFeatures(Features, []);
    end
end         
if (isfield(outer_pars, 'print_every') && ...
	mod(iteration_num, outer_pars.print_every) == 0)
    fprintf('Printing.\n');
	eval(sprintf('print -dpsc2 -append /users-lcv/chaitu/cbp_timit2/figures/%s_features.ps', pars.savename));
end

if (isfield(outer_pars, 'save_every') && ...
	mod(iteration_num, outer_pars.save_every) == 0)
    fprintf('Saving.\n');
    save(sprintf('/users-lcv/chaitu/cbp_timit2/results/%s.mat', pars.savename), ...
         'iteration_num', 'Features', 'outer_pars', 'pars');
end
 
          
% Sliding window length (show the last X iterations).
window_len = 1e3;

% Figure 2 : Plot residual norm distribution
sfigure(2); clf;
PlotResidualStats(data, reconstructed_data, pars.noise_sigma, font_size);
title('Residual distribution');
           
% Figure 3 : Plot magnitude statistics
sfigure(3); clf;
PlotMagnitudeStats(Magnitudes);
    

% Figure 4 : Plot the progress of the L2 term, sparsity term, and
% overall objective.
sfigure(4); clf;
subplot(3, 1, 1);
PlotSparsityProgress(iteration_num, SparsityProgress, window_len,...
                     font_size, line_width, marker_size)

subplot(3, 1, 2);
PlotL2Progress(iteration_num, ResidualNormProgress, window_len,...
               font_size, line_width, marker_size);
subplot(3, 1, 3);
idx = max(1, iteration_num - window_len) : iteration_num; % iterations to plot
plot(idx,ResidualNormProgress(idx) + SparsityProgress(idx), '.-', ...
     'LineWidth', line_width, 'MarkerSize', marker_size);
if (iteration_num > 1)
    xlim([idx(1) idx(end)]);
end
set(gca,'FontSize', font_size);
xlabel('iteration number');
ylabel('objective value');                                  
DisplayObjectiveValueInTitle(iteration_num, ...
                             ResidualNormProgress, SparsityProgress);

% Figure 5 : Plot a randomly chosen sample reconstructed/event train 
sfigure(5); clf;
sample_num = randsample(length(data), 1);
visualize_soln(TransformParams{sample_num}, Magnitudes{sample_num}, ...
               data{sample_num}, reconstructed_data{sample_num}, ...
               0, pars);

           
pause(0.25);