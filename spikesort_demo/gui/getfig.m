function hfig = getfig(hobj)
% GETFIG     Return the handle for the ancestor figure that contains the object referred to in HOBJ
%
% usage:  hfig = getfig(hobj)
%
% arguments:  hobj - Handle for the object whose ancestor figure we want to find
%
% 2010-01 phli
%

% Iteratively ascend handle hierarchy until we get to a figure
hcurr = hobj;
while ~strcmp(get(hcurr, 'Type'), 'figure')
    hcurr = get(hcurr, 'Parent');
    
    if isempty(hcurr)
        error('GETFIG:nofigureancestor', 'No figure found in ancestor hierarchy for given handle');
    end
end

hfig = hcurr;