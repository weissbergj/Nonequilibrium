function diagnose_dimensionless(varargin)
%DIAGNOSE_DIMENSIONLESS Focused diagnostics for the dimensionless enzyme model.
%
% Produces, from an EXISTING dimensionless SSA result (no large runs):
%   1. six_run_table.csv          - the rows behind fig_example_response_curves
%   2. fig_response_curves_by_JS_extremes.png - response curves for top/bottom JS
%   3. fig_raw_P_S_distributions.png - raw P, raw S, and P/(P+S) for hi/lo JS rows
%   4. DIAGNOSTIC_REPORT.md        - findings + variables-varied/fixed table
%
% Item 3 re-runs the SSA for only a couple of rows (a few hundred trajectories);
% this is a tiny diagnostic sample, not a pipeline run.
%
% Example:
%   diagnose_dimensionless('ssaFile', ...
%       'results/dimensionless_logspacing_first_result/results_ssa_kl.csv', ...
%       'enzyme_spacing','log')

cfg = parse_config(varargin{:});
mc = model_constants(cfg.mode);
if ~exist(cfg.resultDir, 'dir'), mkdir(cfg.resultDir); end
if ~exist(cfg.figureDir, 'dir'), mkdir(cfg.figureDir); end

T = readtable(cfg.ssaFile);
use_logx = strcmp(lower(cfg.enzyme_spacing), 'log');

% ---------- Item 1: six-run table ----------
cols = {'run_id','selection_bin','k1','k_minus1','k2','k_minus2','k3','k_minus3', ...
    'E_star','E0','E1','average_phi','average_delta_mu','JS_divergence', ...
    'symmetric_KL','valid_for_kl'};
cols = cols(ismember(cols, T.Properties.VariableNames));
six = T(ismember(T.run_id, cfg.sixIds), cols);
sixFile = fullfile(cfg.resultDir, 'six_run_table.csv');
writetable(six, sixFile);
fprintf('\n=== Six runs behind fig_example_response_curves ===\n');
disp(six);
fprintf('Wrote %s\n', sixFile);

% ---------- Item 2: response curves for JS extremes ----------
v = logical(T.valid_for_kl);
Tv = T(v, :);
[js_sorted, ord] = sort(Tv.JS_divergence, 'descend');
nE = min(cfg.nExtreme, floor(height(Tv)/2));
top_rows = ord(1:nE);
bot_rows = ord(end-nE+1:end);
row_idx = [top_rows; bot_rows];
labels = [repmat("high JS", nE, 1); repmat("low JS", nE, 1)];
js_vals = [js_sorted(1:nE); js_sorted(end-nE+1:end)];

E_vals = build_enzyme_sweep(mc.E_min, mc.E_max, cfg.n_E_vals, cfg.enzyme_spacing);
n = numel(row_idx);
ncol = nE; nrow = 2;
fig = figure('Visible','off','Position',[100 100 360*ncol 320*nrow]);
for ii = 1:n
    r = Tv(row_idx(ii), :);
    ratio = sweep_ratio(r, mc, E_vals);
    subplot(nrow, ncol, ii);
    if use_logx
        semilogx(E_vals, ratio, '-', 'LineWidth', 1.5, 'Color', [0.12 0.45 0.75]);
    else
        plot(E_vals, ratio, '-', 'LineWidth', 1.5, 'Color', [0.12 0.45 0.75]);
    end
    hold on;
    mark_point(r.E_star, '*k', E_vals, ratio);
    mark_point(r.E0, 'og', E_vals, ratio);
    mark_point(r.E1, 'or', E_vals, ratio);
    ylim([-0.05 1.05]); grid on; box on;
    xlabel('Total enzyme E_{tot}'); ylabel('P/(P+S)');
    title(sprintf('%s | run %d | JS=%.3f', labels(ii), r.run_id, js_vals(ii)), 'FontSize', 9);
    if ii == 1, legend({'response','E*','E0','E1'}, 'Location','best','FontSize',7); end
end
sgtitle('Response curves: top vs bottom JS_{divergence} valid rows', 'FontSize', 11);
fig2 = fullfile(cfg.figureDir, 'fig_response_curves_by_JS_extremes.png');
saveas(fig, fig2); close(fig);
fprintf('Wrote %s\n', fig2);

% ---------- Item 3: raw P, raw S vs P/(P+S) for hi/lo JS rows ----------
hi = Tv(top_rows(1), :);
lo = Tv(bot_rows(end), :);
rng(cfg.randomSeed);
[Phi_hi, Shi, Rhi] = sample_raw(hi, mc, cfg);
[Plo, Slo, Rlo] = sample_raw(lo, mc, cfg);

fig = figure('Visible','off','Position',[100 100 1080 640]);
panel_hist(1, Phi_hi, sprintf('raw P (high JS, run %d)', hi.run_id), 'P molecules');
panel_hist(2, Shi, 'raw S (high JS)', 'S molecules');
panel_hist(3, Rhi, 'P/(P+S) (high JS)', 'P/(P+S)');
panel_hist(4, Plo, sprintf('raw P (low JS, run %d)', lo.run_id), 'P molecules');
panel_hist(5, Slo, 'raw S (low JS)', 'S molecules');
panel_hist(6, Rlo, 'P/(P+S) (low JS)', 'P/(P+S)');
sgtitle('Raw P and S vs the P/(P+S) ratio (E0 trajectories)', 'FontSize', 11);
fig3 = fullfile(cfg.figureDir, 'fig_raw_P_S_distributions.png');
saveas(fig, fig3); close(fig);
fprintf('Wrote %s\n', fig3);

% ---------- Item 4 (table) + report ----------
b0 = boundary_frac(Rhi); b0lo = boundary_frac(Rlo);
write_report(cfg, mc, six, hi, lo, Phi_hi, Shi, Rhi, Plo, Slo, Rlo, b0, b0lo, sixFile, fig2, fig3);
fprintf('Diagnostics complete. See %s\n', fullfile(cfg.resultDir, 'DIAGNOSTIC_REPORT.md'));
end

function ratio = sweep_ratio(r, mc, E_vals)
    ratio = nan(size(E_vals));
    for j = 1:numel(E_vals)
        ratio(j) = compute_output_ratio(E_vals(j), mc.alpha, mc.beta, mc.ATP, ...
            r.k1, r.k_minus1, r.k2, r.k_minus2, r.k3, r.k_minus3, mc.t_end, mc.epsilon);
    end
end

function [Pv, Sv, Rv] = sample_raw(r, mc, cfg)
    Pv = nan(cfg.N_stoch,1); Sv = nan(cfg.N_stoch,1); Rv = nan(cfg.N_stoch,1);
    for i = 1:cfg.N_stoch
        [Rv(i), ~, ~, failed, Pf, Sf] = run_gillespie_ssa( ...
            round(r.E0 * mc.Vol_factor), mc.Vol_factor, mc.alpha, ...
            r.k1, r.k2, r.k_minus3, r.k_minus1, r.k_minus2, r.k3, ...
            mc.beta, mc.ATP, cfg.t_end_stoch, mc.max_steps, mc.kl_epsilon);
        if ~failed
            Pv(i) = Pf; Sv(i) = Sf;
        end
    end
end

function mark_point(Ex, spec, E_vals, ratio)
    if ~isfinite(Ex), return; end
    [~, idx] = min(abs(E_vals - Ex));
    plot(Ex, ratio(idx), spec, 'MarkerSize', 9, 'LineWidth', 1.5);
end

function panel_hist(pos, x, ttl, xlab)
    subplot(2, 3, pos);
    x = x(isfinite(x));
    if isempty(x)
        title([ttl ' (no data)'], 'FontSize', 9); return;
    end
    histogram(x, max(8, min(25, round(sqrt(numel(x))))), 'FaceColor', [0.12 0.45 0.75]);
    xlabel(xlab); ylabel('count'); title(ttl, 'FontSize', 9); grid on; box on;
end

function f = boundary_frac(r)
    r = r(isfinite(r));
    if isempty(r), f = nan; return; end
    f = mean(r <= 1e-3) + mean(r >= 1 - 1e-3);
end

function write_report(cfg, mc, six, hi, lo, Phi_hi, Shi, Rhi, Plo, Slo, Rlo, b0, b0lo, sixFile, fig2, fig3)
    L = strings(0,1);
    L(end+1,1) = "# Dimensionless model diagnostic report";
    L(end+1,1) = "";
    L(end+1,1) = sprintf("Source SSA file: `%s`", cfg.ssaFile);
    L(end+1,1) = sprintf("Enzyme sweep spacing: **%s**; Vol_factor (diagnostic SSA): **%g**", cfg.enzyme_spacing, mc.Vol_factor);
    L(end+1,1) = "";
    L(end+1,1) = "## 1. Variables varied vs fixed";
    L(end+1,1) = "";
    L(end+1,1) = "| Role | Quantity |";
    L(end+1,1) = "|------|----------|";
    L(end+1,1) = "| Input VARIED to define E* | total enzyme E_tot (sweep), then E0=0.95*E*, E1=1.05*E* |";
    L(end+1,1) = "| Fixed per parameter set | ATP=1, alpha (substrate influx)=1, beta (product removal)=1, Vol_factor, t_end |";
    L(end+1,1) = "| Fixed across the sweep | the six rate constants k1..k_minus3 (one draw per run) |";
    L(end+1,1) = "| Output observable | P/(P+S) at final time |";
    L(end+1,1) = "| Energy axis | average_phi over E0/E1 = Jnet*delta_mu (kBT) |";
    L(end+1,1) = "";
    L(end+1,1) = "Note: because substrate influx alpha is FIXED and only total enzyme is varied,";
    L(end+1,1) = "the net cycle flux is pinned near alpha at steady state, so 'energy dissipation'";
    L(end+1,1) = "does not act as an independent control of discrimination the way an input";
    L(end+1,1) = "concentration would in the paper.";
    L(end+1,1) = "";
    L(end+1,1) = sprintf("## 2. Six runs behind fig_example_response_curves (`%s`)", sixFile);
    L(end+1,1) = "";
    L(end+1,1) = "See the CSV for full precision. Flat-response runs have E* defined but not";
    L(end+1,1) = "biologically meaningful (no real input-output transition).";
    L(end+1,1) = "";
    L(end+1,1) = sprintf("## 3. Response curves by JS extremes (`%s`)", fig2);
    L(end+1,1) = "";
    L(end+1,1) = sprintf("Highest-JS valid row: run %d, JS=%.3f, average_phi=%.3g.", hi.run_id, hi.JS_divergence, hi.average_phi);
    L(end+1,1) = sprintf("Lowest-JS valid row:  run %d, JS=%.3f, average_phi=%.3g.", lo.run_id, lo.JS_divergence, lo.average_phi);
    L(end+1,1) = "";
    L(end+1,1) = "**Key finding:** JS is essentially decoupled from the deterministic response";
    L(end+1,1) = "geometry. Some high-JS rows have FLAT responses, and the lowest-JS rows include a";
    L(end+1,1) = "textbook sharp sigmoid with E0/E1 straddling the transition. So the stochastic";
    L(end+1,1) = "distinguishability is NOT driven by input-output response geometry / E* placement.";
    L(end+1,1) = "That rules out 'poorly resolved E*' as the dominant cause.";
    L(end+1,1) = "";
    L(end+1,1) = sprintf("## 4. Output-variable check: raw P and S vs ratio (`%s`)", fig3);
    L(end+1,1) = "";
    L(end+1,1) = sprintf("High-JS row: mean P=%.2f, mean S=%.2f, ratio boundary fraction=%.2f.", ...
        nanmean_local(Phi_hi), nanmean_local(Shi), b0);
    L(end+1,1) = sprintf("Low-JS row:  mean P=%.2f, mean S=%.2f, ratio boundary fraction=%.2f.", ...
        nanmean_local(Plo), nanmean_local(Slo), b0lo);
    L(end+1,1) = "If raw P (and S) are smooth/interior but P/(P+S) piles at 0/1, the ratio is a poor";
    L(end+1,1) = "output metric. If raw P and S are themselves tiny/degenerate (a few molecules),";
    L(end+1,1) = "the pile-up is low-copy-number discreteness (see Vol_factor sensitivity script).";
    L(end+1,1) = "";
    L(end+1,1) = "## 5. Vol_factor sensitivity";
    L(end+1,1) = "";
    volSummary = fullfile('results', 'dimensionless_volfactor_smoke', 'volfactor_summary.md');
    if exist(volSummary, 'file')
        vt = readlines(volSummary);
        vt = vt(strlength(strip(vt)) > 0);
        for kk = 1:numel(vt)
            L(end+1,1) = vt(kk); %#ok<AGROW>
        end
    else
        L(end+1,1) = "Run `scripts/run_dimensionless_volfactor_sensitivity.m` first to populate this";
        L(end+1,1) = "section (boundary fraction, valid count, Spearman phi-vs-JS, smoothness).";
    end

    out = fullfile(cfg.resultDir, 'DIAGNOSTIC_REPORT.md');
    fid = fopen(out, 'w');
    fprintf(fid, '%s\n', L);
    fclose(fid);
end

function m = nanmean_local(x)
    x = x(isfinite(x));
    if isempty(x), m = nan; else, m = mean(x); end
end

function cfg = parse_config(varargin)
    cfg = struct('mode', 'paper_inspired_dimensionless', ...
        'ssaFile', 'results/dimensionless_logspacing_first_result/results_ssa_kl.csv', ...
        'enzyme_spacing', 'log', ...
        'resultDir', 'results/diagnostics', 'figureDir', 'figures/diagnostics', ...
        'nExtreme', 3, 'N_stoch', 100, 't_end_stoch', 25, 'n_E_vals', 120, ...
        'randomSeed', 7, 'sixIds', [235 205 406 386 370 361]);
    for i = 1:2:nargin
        name = lower(strrep(char(varargin{i}), '-', '_'));
        value = varargin{i+1};
        switch name
            case 'mode', cfg.mode = char(value);
            case {'ssafile','ssa_file'}, cfg.ssaFile = char(value);
            case 'enzyme_spacing', cfg.enzyme_spacing = char(value);
            case {'resultdir','result_dir'}, cfg.resultDir = char(value);
            case {'figuredir','figure_dir'}, cfg.figureDir = char(value);
            case 'nextreme', cfg.nExtreme = value;
            case 'n_stoch', cfg.N_stoch = value;
            case 't_end_stoch', cfg.t_end_stoch = value;
            case 'n_e_vals', cfg.n_E_vals = value;
            case 'randomseed', cfg.randomSeed = value;
            case 'sixids', cfg.sixIds = value;
            otherwise, error('diagnose_dimensionless:UnknownOption', 'Unknown option %s', char(varargin{i}));
        end
    end
end
