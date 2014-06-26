function fig = gen_slider_plot(varargin)
%GEN_SLIDER_PLOT     Basic set up for a plot with a slider; usually called from another more specific function
% Copied from SNL-E matlab-standard and stripped down

opts = inputParser();
opts.addParamValue('handle', []);
opts.addParamValue('data', []);
opts.addParamValue('data_length_func', @length);
opts.addParamValue('slider_func',      @slider_func);
opts.addParamValue('get_data_func',    @get_data_func);
opts.addParamValue('plot_func',        @plot_func);
opts.parse(varargin{:});
opts = opts.Results;

if isempty(opts.handle)
    opts.handle = figure();
end
figure(opts.handle);
fig = getfig(opts.handle);

% Add data to figure
setappdata(fig, 'data', opts.data);


% Create an API to simplify calling GUI methods from outside
api = getappdata(fig, 'api');
if isempty(api)
    api = struct();
end
api.data_length = opts.data_length_func;
api.slider_plot = opts.slider_func;
api.get_data    = opts.get_data_func;
api.plot        = opts.plot_func;
setappdata(fig, 'api', api);


% Add image panel to fig
% impanel = uipanel('Units', 'normalized', 'Position', [0, 0, 1, 1], 'Parent', opts.handle, 'Tag', 'impanel');


% For slider, switch to Metal Look and Feel if running on Mac OS > 10.6
nativelaf = javax.swing.UIManager.getLookAndFeel();
if aftersnowleopard()
    javax.swing.UIManager.setLookAndFeel('javax.swing.plaf.metal.MetalLookAndFeel');
end

% Add slider to panel
start_index = 1;
index_min = 0.99; 
index_max = api.data_length(opts.data);
slider_step = 1 / (index_max - index_min);
slider_steps = min([slider_step 1]) * [1 5];
shandle = uicontrol(fig, ...
    'Style'     , 'slider',                        ...
    'Min'       , index_min,                       ...
    'Max'       , index_max,                       ...
    'Units'     , 'normalized',                    ...
    'Position'  , [0.01, 0.025, 0.95, 0.01],       ...
    'Value'     , start_index,                     ...
    'SliderStep', slider_steps,                    ...
    'CallBack'  , api.slider_plot,                ...
    'Tag'       , 'slider'                   ...
);
drawnow;
javax.swing.UIManager.setLookAndFeel(nativelaf);
setappdata(fig, 'slider', shandle);

% Control with keyboard arrow keys
if haveipt()
    add_slider_arrowkeys(fig, shandle);
end


% Set up figure
set(opts.handle, 'Toolbar', 'figure')

% Set up axes in panel
% ax = axes('Parent', fig, 'Units', 'normalized', 'Position', [0, 0.025, 0.95, 0.9]);
ax = axes();
setappdata(fig, 'slider_axes', ax);

% Initial plot call
api.slider_plot(shandle, []);

% Now set axes to replace children
set(ax, 'NextPlot', 'replaceChildren');


if nargout < 1, clear fig; end


function slider_func(slider, ev) %#ok<INUSD>
f = getfig(slider);
api = getappdata(f, 'api');

% Determine what value the slider is set to
sval = round(get(slider, 'Value'));
pval = round(getappdata(slider, 'PreviousValue'));
if sval == pval, return; end

% Set axes
ax = getappdata(f, 'slider_axes');
axes(ax);

% Plot image data from stack
data = getappdata(f, 'data');
d = api.get_data(data, sval);
i = api.plot(d, sval, data, f);

% Set handle
setappdata(f, 'imhandle', i);

% Useful for higher level guis that add callbacks before this one, e.g. stack_point_picker
setappdata(slider, 'PreviousValue', sval);


function h = plot_func(d, sval, data, f)
h = plot(d);
title(num2str(sval));

function d = get_data_func(data, sval)
d = data{sval};