function repoRoot = setup_repo_paths()
%SETUP_REPO_PATHS Add pipeline paths when called from scripts/ entry points.

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
cd(repoRoot);

addpath(genpath(fullfile(repoRoot, 'src', 'core')));
addpath(genpath(fullfile(repoRoot, 'src', 'pipeline')));
addpath(genpath(fullfile(repoRoot, 'src', 'diagnostics')));

fprintf('Nonequilibrium repo root: %s\n', repoRoot);

end
