function add_slider_arrowkeys(f, slider)
% Copied from SNL-E matlab-standard
% phli
%

iptaddcallback(f, 'WindowKeyPressFcn', @(fig, ev) slidearrowkeys(ev, slider));

function slidearrowkeys(ev, slider)
switch ev.Key
    case {'leftarrow',  'uparrow'}
        if any(cellfun(@(s) strcmp(s, 'shift'), ev.Modifier)), incr = -5;
        else                                                   incr = -1; end
    case {'rightarrow', 'downarrow'}
        if any(cellfun(@(s) strcmp(s, 'shift'), ev.Modifier)), incr = 5;
        else                                                   incr = 1; end
    otherwise
        return
end

sval = round(get(slider, 'Value'));

% Check boundaries
newval = sval + incr;
newval = max(newval, get(slider, 'Min'));
newval = min(newval, get(slider, 'Max'));

set(slider, 'Value', newval);
cb = get(slider, 'Callback');
cb(slider, []);