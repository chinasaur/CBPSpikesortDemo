function PlotPartition(signal, IDX, pars, plot_len, fig1)
    if (~exist('fig1', 'var'))
        fig1 = figure;
    end
    if (~exist('plot_len', 'var'))
        plot_len = 5e2;
    end
    middle_idx = ceil(size(signal, 2) / 2) + (-plot_len : plot_len);
    colors = hsv(5);
    figure(fig1); clf; hold on;
    tax = (0:size(signal,2)-1) .* 1/16;
    tax = tax - tax(middle_idx(1));
    signal_rms = sqrt(sum(signal .^ 2, 1));    
    plot(tax, signal_rms, 'k-');        
    k = 1;
    
    1;
    for i = 1 : length(IDX)
        if (IDX{i}(end) < middle_idx(1) || IDX{i}(1) > middle_idx(end))
            continue;
        end
        idx = IDX{i};
        plot(tax(idx), reshape(signal_rms(idx), [], 1), ...
             'Color', colors(mod(k - 1, size(colors, 1)) + 1, :));  
        k = k + 1;
    end
    xl = [tax(middle_idx(1)), tax(middle_idx(end))];
    xlim(xl);    
    plot(xl, pars.threshold .* [1 1], 'k');   
    font_size = 16;
    set(gca, 'FontSize', font_size);
    xlabel('time (ms)');
    ylabel('trace');    
end
