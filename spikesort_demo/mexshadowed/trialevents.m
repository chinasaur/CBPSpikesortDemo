function [dts inds] = trialevents(events, triggers, t0, t1)

warning('TRIALEVENTS:no_mex', 'Running slow uncompiled version; please compile mex/trialevents.c');

n = numel(triggers);
events = events(:);
dts = cell(n,1);
inds = cell(n,1);
for i = 1:n
    dt = events - triggers(i);
    dt = dt(dt >= t0 & dt <= t1);
    dts{i} = dt;
    inds{i} = i * ones(length(dt),1);
end
dts = cell2mat(dts);
inds = cell2mat(inds);