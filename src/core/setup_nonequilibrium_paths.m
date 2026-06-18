function repoRoot = setup_nonequilibrium_paths()
%SETUP_NONEQUILIBRIUM_PATHS Add src folders to the MATLAB path from repo root.
%
%   repoRoot = setup_nonequilibrium_paths()
%
% Call from scripts/run_*.m after changing to the repository root.

thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(thisFile)));

addpath(genpath(fullfile(repoRoot, 'src', 'core')));
addpath(genpath(fullfile(repoRoot, 'src', 'pipeline')));
addpath(genpath(fullfile(repoRoot, 'src', 'diagnostics')));

end
