function spikesort_demo_setup(ssdpath)

addpath(genpath(ssdpath))
%addpath(fullfile(ssdpath, '../example_data/'));
javaaddpath(fullfile(ssdpath, 'cbp/utils/parforProgress/'));

% Check that MEX files are compiled (the following lines will print
% warnings if not).
[~, ~] = greedymatchtimes([], [], [], []);
[~, ~] = trialevents([], [], [], []);
ecos(1, sparse(1), 1, struct('l', 1, 'q', []), struct('verbose', 0))

% Setup matlabpool for parfor if available and no pool open already
if(0)
if ((exist('matlabpool')==2) && (matlabpool('size') == 0))
  try
    matlabpool open
  catch me
    warning('Failed to open parallel sessions using matlabpool:\n  %s\n',...
        me.message);
  end
end
end

if (exist('parpool')==2)
  try
    if (isempty(gcp('nocreate')))
      parpool
    end
  catch me
    warning('Failed to open parallel pool using parpool:\n  %s\n',...
        me.message);
  end
end
