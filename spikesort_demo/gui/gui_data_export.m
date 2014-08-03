function submenu = gui_data_export(handle, dataname, menu_name_or_handle)
% GUI_DATA_EXPORT    Set up a menu item for saving data
% usage: submenu = gui_data_export(handle, dataname, menu_name_or_handle)
%
% Sets up a menu item in the figure window for given handle that exports 
% the appdata corresponding to the given dataname to the base workspace.
% The item can either be added to an existing top level menu, in which case
% the top menu handle is passed, or else the item can be added into a new
% top level menu, in which case the menu name is passed as a string.
%
% Returns the handle to the created submenu.
%
% 2010-05 phli
%

fig = getfig(handle);

% Get existing top menu or create
if ischar(menu_name_or_handle)
    menu_name = menu_name_or_handle;
    menu = uimenu(fig, 'Label', menu_name, 'Tag', menu_name);
else
    menu = menu_name_or_handle;
end

% Add submenu item to save data
submenu = uimenu(menu, 'Label', ['Export ' dataname ' to workspace'], 'Tag', ['save_' dataname '_menu'], 'Callback', {@export_data, dataname});

% Add export_data to api
api = getappdata(fig, 'api');
if isempty(api)
    api = struct();
end
api.export_data = @export_data;
api.(['export_' dataname]) = get(submenu, 'Callback');
setappdata(fig, 'api', api);


function export_data(menu, ev, dataname) %#ok<INUSL>
fig = getfig(menu);
data = getappdata(fig, dataname);
assignin('base', dataname, data);