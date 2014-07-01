function [f thresh] = AmplitudeThresholdGUI(spikeamps, spiketimes, varargin)
opts = inputParser();
opts.addParamValue('f', []);
opts.addParamValue('thresh', []);
opts.addParamValue('ampbins', 1e2);
%opts.addParamValue('fontsize', 16);
opts.addParamValue('kdepoints', 32);
opts.addParamValue('kderange', [0.3 1.1]);
opts.addParamValue('kdewidth', 5);
opts.addParamValue('dt', 1);
opts.addParamValue('true_sp', {});
opts.addParamValue('location_slack', 30); % Allowed mismatch in timing for comparison to ground truth
opts.parse(varargin{:});
opts = opts.Results;

% Setup figure
f = opts.f;
if isempty(f), f = figure(); end

n = length(spikeamps);
setappdata(f, 'spikeamps', spikeamps);

% Calculate KDE thresholds if none given
thresh = opts.thresh;
if isempty(thresh)
    for i = 1:n
        thresh(i) = ComputeKDEThreshold(spikeamps{i}, opts.kdepoints, opts.kderange, opts.kdewidth);
    end
end

% Store initial thresh value
guiout = getappdata(0, 'GUIOutput');
guiout{f}.thresh = thresh;
setappdata(0, 'GUIOutput', guiout);

% Modify spiketimes by dt
spiketimes = cellfun(@(st) st.*opts.dt, spiketimes,   'UniformOutput', false);
true_sp    = cellfun(@(st) st.*opts.dt, opts.true_sp, 'UniformOutput', false);
slack = opts.location_slack*opts.dt;
setappdata(f, 'spiketimes', spiketimes);
setappdata(f, 'true_sp', true_sp);
setappdata(f, 'location_slack', slack);

% Store initial thresholding
threshspiketimes = cell(size(spiketimes));
for i = 1:n
    threshspiketimes{i} = spiketimes{i}(spikeamps{i} > thresh(i));
end
setappdata(f, 'threshspiketimes', threshspiketimes);

v = ver();
haveipt = any(strcmp('Image Processing Toolbox', {v.Name}));
figure(f);
for i = 1:n
    subplot(n+1, n, i);
    
    % Plot spike amplitude histogram
    hist(spikeamps{i}, opts.ampbins);
%    set(gca, 'FontSize', opts.fontsize);
    title(sprintf('Amplitudes, cell %d', i));
    xlim([0 max([spikeamps{i}(:); 1.5])]);
    
    % Plot threshold as vertical red lines.
    hold on;    
    yl = get(gca, 'YLim');
    if haveipt
        xl = get(gca, 'XLim');
        cnstrfcn = makeConstrainToRectFcn('imline', xl, yl);
        lh = imline(gca, thresh(i)*[1 1], yl, 'PositionConstraintFcn', cnstrfcn);
        lh.setColor('r');
        lch = get(lh, 'Children');
        set(lch(1:2), 'HitTest', 'off');
        lh.addNewPositionCallback(@(pos) updateThresh(pos(1), i, f));
    else
        plot(thresh(i) * [1 1], yl, 'r-', 'LineWidth', 2);
    end
    ylim(yl);
end


% Plot initial ACorr/XCorrs
for i = 1:n
    plotACorr(threshspiketimes, i);
    for j = (i+1):n
        plotXCorr(threshspiketimes, i, j);
    end
end


% Report on performance relative to ground truth if available
showGroundTruthEval(threshspiketimes, f);

if nargout < 1, clear f; end
if nargout > 1
    waitfor(f);
    guiout = getappdata(0, 'GUIOutput');
    thresh = guiout{f}.thresh;
end


function showGroundTruthEval(spiketimes, f)
true_sp = getappdata(f, 'true_sp');
if isempty(true_sp), return; end
slack = getappdata(f, 'location_slack');

% Evaluate CBP sorting
[total_misses, total_false_positives] = ...
    evaluate_sorting(spiketimes, true_sp, slack);

% Display on fig
n = length(spiketimes);
for i = 1:length(true_sp)
    if isempty(true_sp{i}), continue; end
    subplot(n+1, n, i);
    xlabel(sprintf('misses: %d fps: %d', total_misses(i), total_false_positives(i)));
end


function updateThresh(newthresh, i, f)
guiout = getappdata(0, 'GUIOutput');
threshsts = getappdata(f, 'threshspiketimes');
sts       = getappdata(f, 'spiketimes');
amps      = getappdata(f, 'spikeamps');

% Calculate new threshed sts
threshsts{i} = sts{i}(amps{i} > newthresh);

% Plot
plotACorr(threshsts, i);
n = length(threshsts);
for j = 1:n
    if j == i, continue; end
    plotXCorr(threshsts, i, j);
end

showGroundTruthEval(threshsts, f);

% Save new threshes and threshed spiketimes
setappdata(f, 'threshspiketimes', threshsts);
guiout{f}.thresh(i) = newthresh;
setappdata(0, 'GUIOutput', guiout);


function plotACorr(spiketimes, i)
n = length(spiketimes);
subplot(n+1, n, sub2ind([n n+1], i, n+1));
cla;
psthacorr(spiketimes{i})

function plotXCorr(spiketimes, i, j)
if j < i
    tmp = i;
    i = j;
    j = tmp;
end

n = length(spiketimes);
subplot(n+1, n, sub2ind([n n+1], j, i+1));
cla;
psthxcorr(spiketimes{i}, spiketimes{j})