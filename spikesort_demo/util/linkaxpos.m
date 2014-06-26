function linkaxpos(varargin)
% Copied from plotyy

if length(varargin) == 1, ax = varargin{1};
else                      ax = cell2mat(varargin); end

if ~feature('HGUsingMATLABClasses')
    listenerFun = @handle.listener;
    postSetStr = 'PropertyPostSet';
else
    listenerFun = @event.proplistener;
    postSetStr = 'PostSet';
end
hList(1) = listenerFun(handle(ax(1)),findprop(handle(ax(1)),'Position'),...
    postSetStr,@(obj,evd)(localUpdatePosition(obj,evd,ax(1),ax(2))));
hList(2) = listenerFun(handle(ax(2)),findprop(handle(ax(2)),'Position'),...
    postSetStr,@(obj,evd)(localUpdatePosition(obj,evd,ax(2),ax(1))));
setappdata(ax(1),'graphicsPlotyyPositionListener',hList);


% Copied from plotyy
% Keep the positions of two axes in sync:
function localUpdatePosition(~,~,axSource,axDest)
newPos = get(axSource,'Position');
hFig = ancestor(axSource,'Figure');
newDestPos = hgconvertunits(hFig,newPos,get(axSource,'Units'),get(axDest,'Units'),get(axSource,'Parent'));
set(axDest,'Position',newDestPos);