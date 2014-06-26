function isnip = FindSnippets(findsnips, snipcenters, snipwidths)
% FindSnippets
% usage: isnip = FindSnippets(findsnip, snipcenters, snipwidths)
%
% Get indices of snippets containing time points FINDSNIPS.  Returns zero
% for unfound indices.
%
% As shorthand, can also call as:
%        isnip = FindSnippets(findsnip, snipcenters, snippets)
% in which case SNIPWIDTHS are calculated automatically
%

if iscell(snipwidths)
    snipwidths = cellfun(@(S) size(S,1), snipwidths);
end

isnip = zeros(size(findsnips));
for i = 1:numel(findsnips)
    [vmin imin] = min(abs(snipcenters - findsnips(i)));
    if vmin <= 0.5 + snipwidths(imin) / 2
        isnip(i) = imin;
    end
end