%RUN_DIMENSIONLESS_VOLFACTOR_SENSITIVITY Tiny opt-in Vol_factor sweep.
%
% From repo root:
%   matlab -batch "run('scripts/run_dimensionless_volfactor_sensitivity.m')"
%
% Question: is the boundary pile-up in P/(P+S) caused by low molecule counts
% (Vol_factor=1), i.e. low-copy-number discreteness, rather than true model
% behavior? Runs ONE small log-spaced screen+selection, then re-runs the SSA at
% Vol_factor = 1, 10, 100 on the SAME selected rows and compares:
%   - median boundary fraction (fraction_at_0 + fraction_at_1)
%   - valid_for_kl count
%   - Spearman(average_phi, JS_divergence)
%   - median output variance (a smoothness proxy)
%
% This is a deliberately small diagnostic, not a pipeline run.

setup_repo_paths();
tStart = tic;

mode = 'paper_inspired_dimensionless';
resultDir = fullfile('results', 'dimensionless_volfactor_smoke');
figureDir = fullfile('figures', 'dimensionless_volfactor_smoke');
if ~exist(resultDir, 'dir'), mkdir(resultDir); end
if ~exist(figureDir, 'dir'), mkdir(figureDir); end

screenFile = fullfile(resultDir, 'results_deterministic_screen.csv');
selectedFile = fullfile(resultDir, 'selected_parameter_sets.csv');

fprintf('\n=== Vol_factor sensitivity (log-spaced, small) ===\n');

% Deliberately small: bounded SSA cost via a low max_steps cap. Higher Vol_factor
% means many more SSA events, so without the cap a Vol_factor=100 run becomes a
% large run. Rows that hit the cap are reported (they are NOT valid_for_kl).
N_select = 8;
N_stoch = 80;        % must exceed 5*nBins so rows are not flagged thin
max_steps_cap = 2e5;

run_deterministic_screen('mode', mode, 'enzyme_spacing', 'log', ...
    'N_runs', 60, 'n_E_vals', 60, 'useParallel', false, ...
    'outputFile', screenFile, 'randomSeed', 1);

select_parameter_sets('inputFile', screenFile, 'N_select', N_select, ...
    'oversampleFactor', 1.0, 'outputFile', selectedFile, ...
    'randomSeed', 2, 'stratifyBy', 'average_phi');

Vfs = [1 10 100];
rows = strings(0,1);
rows(end+1,1) = "# Vol_factor sensitivity (dimensionless, log-spaced)";
rows(end+1,1) = "";
rows(end+1,1) = sprintf("Bounded test: N_select=%d, N_stoch=%d, max_steps cap=%g.", N_select, N_stoch, max_steps_cap);
rows(end+1,1) = "";
rows(end+1,1) = "| Vol_factor | valid_for_kl | max_steps_hit | median boundary frac | Spearman phi-vs-JS | median output_var |";
rows(end+1,1) = "|-----------:|-------------:|--------------:|---------------------:|-------------------:|------------------:|";

for vi = 1:numel(Vfs)
    V = Vfs(vi);
    ssaFile = fullfile(resultDir, sprintf('results_ssa_kl_vol%d.csv', V));
    run_ssa_kl_on_selected('mode', mode, 'inputFile', selectedFile, ...
        'Vol_factor', V, 'max_steps', max_steps_cap, 'N_stoch', N_stoch, ...
        't_end_stoch', 25, 'nBins', 10, 'useParallel', false, 'batchSize', 4, ...
        'outputFile', ssaFile, 'randomSeed', 3, 'saveExampleDistributions', false);

    T = readtable(ssaFile);
    v = logical(T.valid_for_kl);
    Tv = T(v, :);
    nValid = height(Tv);
    nMaxsteps = sum(T.max_steps_hits_E0 > 0 | T.max_steps_hits_E1 > 0);
    if nValid > 0
        bf = T.fraction_at_0_E0(v) + T.fraction_at_1_E0(v);
        medBoundary = median(bf(isfinite(bf)));
        medVar = median(T.output_var_E0(v), 'omitnan');
        rho = spearman_local(Tv.average_phi, Tv.JS_divergence);
    else
        medBoundary = nan; medVar = nan; rho = nan;
    end
    rows(end+1,1) = sprintf("| %d | %d | %d | %.2f | %.3f | %.4f |", V, nValid, nMaxsteps, medBoundary, rho, medVar); %#ok<SAGROW>
    fprintf('Vol_factor=%d: valid=%d, max_steps_hit=%d, medBoundary=%.2f, rho=%.3f, medVar=%.4f\n', ...
        V, nValid, nMaxsteps, medBoundary, rho, medVar);
end

rows(end+1,1) = "";
rows(end+1,1) = "Interpretation:";
rows(end+1,1) = "- A sharp drop in median boundary fraction as Vol_factor increases means the 0/1";
rows(end+1,1) = "  pile-up at Vol_factor=1 is a LOW-COPY-NUMBER (discreteness) artifact, not true";
rows(end+1,1) = "  model behavior. (median output_var also shrinks because the ratio concentrates";
rows(end+1,1) = "  near its deterministic value as molecule counts grow.)";
rows(end+1,1) = "- Spearman phi-vs-JS is computed on only the few valid rows here, so its sign/size";
rows(end+1,1) = "  is NOT statistically reliable; treat it as a hint, not evidence of a trend.";
rows(end+1,1) = "  A real phi-vs-JS trend must be confirmed on a full run at the larger Vol_factor.";

summaryFile = fullfile(resultDir, 'volfactor_summary.md');
fid = fopen(summaryFile, 'w');
fprintf(fid, '%s\n', rows);
fclose(fid);

fprintf('\nVol_factor sensitivity finished in %.1f minutes. Wrote %s\n', toc(tStart)/60, summaryFile);

% ---- local functions ----
function rho = spearman_local(x, y)
    x = x(:); y = y(:);
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    if numel(x) < 3, rho = nan; return; end
    rx = rank_with_ties(x); ry = rank_with_ties(y);
    rx = rx - mean(rx); ry = ry - mean(ry);
    denom = sqrt(sum(rx.^2) * sum(ry.^2));
    if denom == 0, rho = nan; else, rho = sum(rx .* ry) / denom; end
end

function r = rank_with_ties(x)
    [xs, order] = sort(x);
    n = numel(x);
    r = zeros(n, 1);
    r(order) = 1:n;
    i = 1;
    while i <= n
        j = i;
        while j < n && xs(j+1) == xs(i), j = j + 1; end
        if j > i, r(order(i:j)) = mean(i:j); end
        i = j + 1;
    end
end
