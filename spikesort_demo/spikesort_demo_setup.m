function spikesort_demo_setup(ssdpath)
addpath(genpath(ssdpath))
%addpath(fullfile(ssdpath, '../example_data/'));
javaaddpath(fullfile(ssdpath, 'cbp/utils/parforProgress/'));

% Setup matlabpool for parfor if available and no pool open already
if exist('matlabpool', 'file') && matlabpool('size') == 0
  try
    matlabpool open
  catch me
    warning('Failed to open parallel sessions using matlabpool:\n  %s\n',...
        me.message);
  end
end

% Check installation
% The following lines will print warnings if there are mex files that need 
% to be compiled.
%
[~, ~] = greedymatchtimes([], [], [], []);
[~, ~] = trialevents([], [], [], []);
ecos(1, sparse(1), 1, struct('l', 1, 'q', []), struct('verbose', 0))