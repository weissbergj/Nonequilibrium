function plot_response_curves(varargin)
%PLOT_RESPONSE_CURVES Diagnostic: deterministic P/(P+S) vs total enzyme.
%
% Re-solves the deterministic enzyme sweep for a few representative parameter
% sets and plots the output ratio against total enzyme, marking E*, E0 and E1.
% This is a sanity check that E* sits on a monotonic/sigmoidal region of the
% response curve rather than on a flat or arbitrary part.
%
% Example:
%   plot_response_curves('mode','paper_inspired_dimensionless', ...
%       'selectedFile','results/dimensionless_first_result/selected_parameter_sets.csv', ...
%       'outputFile','figures/dimensionless_first_result/fig_example_response_curves.png')

cfg = parse_config(varargin{:});
mc = model_constants(cfg.mode);

if exist(cfg.selectedFile, 'file')
    T = readtable(cfg.selectedFile);
elseif exist(cfg.screenFile, 'file')
    T = readtable(cfg.screenFile);
else
    warning('plot_response_curves:NoInput', 'No selected/screen file found; skipping response curves.');
    return
end

rows = pick_rows(T, cfg.nCurves);
if isempty(rows)
    warning('plot_response_curves:NoRows', 'No usable rows for response curves.');
    return
end

spacing = resolve_spacing(cfg.enzyme_spacing, T);
use_logx = strcmp(spacing, 'log');
E_vals = build_enzyme_sweep(mc.E_min, mc.E_max, cfg.n_E_vals, spacing);
out_dir = fileparts(cfg.outputFile);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

n = numel(rows);
ncol = min(3, n);
nrow = ceil(n / ncol);
fig = figure('Visible', 'off', 'Position', [100 100 360*ncol 320*nrow]);

for ii = 1:n
    r = T(rows(ii), :);
    k = [r.k1, r.k_minus1, r.k2, r.k_minus2, r.k3, r.k_minus3];
    ratio = nan(size(E_vals));
    for j = 1:numel(E_vals)
        ratio(j) = compute_output_ratio(E_vals(j), mc.alpha, mc.beta, mc.ATP, ...
            k(1), k(2), k(3), k(4), k(5), k(6), mc.t_end, mc.epsilon);
    end

    subplot(nrow, ncol, ii);
    if use_logx
        semilogx(E_vals, ratio, '-', 'LineWidth', 1.5, 'Color', [0.12 0.45 0.75]);
    else
        plot(E_vals, ratio, '-', 'LineWidth', 1.5, 'Color', [0.12 0.45 0.75]);
    end
    hold on;
    mark_point(r, 'E_star', '*k', E_vals, ratio);
    mark_point(r, 'E0', 'og', E_vals, ratio);
    mark_point(r, 'E1', 'or', E_vals, ratio);
    xlabel('Total enzyme E_{tot}'); ylabel('P/(P+S)');
    title(curve_title(T, rows(ii)), 'Interpreter', 'none', 'FontSize', 9);
    ylim([-0.05 1.05]); grid on; box on;
    if ii == 1
        legend({'response','E*','E0','E1'}, 'Location', 'best', 'FontSize', 7);
    end
end

sgtitle(sprintf('Deterministic response curves (%s, %s spacing): P/(P+S) vs total enzyme', ...
    mc.mode, spacing), 'Interpreter', 'none', 'FontSize', 11);
annotation(fig, 'textbox', [0.02 0.0 0.96 0.045], ...
    'String', 'Note: some selected parameter sets do not exhibit smooth midpoint-like response curves.', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 9, 'FontAngle', 'italic', 'Color', [0.25 0.25 0.25]);
saveas(fig, cfg.outputFile);
close(fig);
fprintf('Wrote response-curve diagnostic to %s (%d curves)\n', cfg.outputFile, n);
end

function mark_point(r, col, spec, E_vals, ratio)
    if ~ismember(col, r.Properties.VariableNames), return; end
    Ex = r.(col);
    if ~isfinite(Ex), return; end
    [~, idx] = min(abs(E_vals - Ex));
    plot(Ex, ratio(idx), spec, 'MarkerSize', 9, 'LineWidth', 1.5);
end

function spacing = resolve_spacing(forced, T)
    if ~isempty(forced)
        spacing = lower(char(forced));
        return
    end
    spacing = 'linear';
    if ismember('enzyme_spacing_setting', T.Properties.VariableNames) && height(T) > 0
        spacing = lower(char(string(T.enzyme_spacing_setting(1))));
    end
end

function ttl = curve_title(T, idx)
    bin = '';
    if ismember('selection_bin', T.Properties.VariableNames)
        bin = char(string(T.selection_bin(idx)));
    end
    rid = idx;
    if ismember('run_id', T.Properties.VariableNames)
        rid = T.run_id(idx);
    end
    if isempty(bin)
        ttl = sprintf('run %d', rid);
    else
        ttl = sprintf('run %d (%s)', rid, bin);
    end
end

function rows = pick_rows(T, nCurves)
    n = height(T);
    if n == 0, rows = []; return; end
    if ismember('selection_bin', T.Properties.VariableNames)
        bins = string(T.selection_bin);
        order = ["low", "medium", "high"];
        rows = [];
        for L = 1:numel(order)
            idx = find(startsWith(bins, order(L)), 1);
            if ~isempty(idx), rows(end+1) = idx; end %#ok<AGROW>
        end
        % Top up to nCurves with evenly spaced extra rows.
        extra = setdiff(round(linspace(1, n, nCurves)), rows, 'stable');
        rows = [rows, reshape(extra, 1, [])];
    else
        rows = unique(round(linspace(1, n, nCurves)), 'stable');
    end
    rows = rows(1:min(nCurves, numel(rows)));
end

function cfg = parse_config(varargin)
    cfg = struct('mode', 'legacy', ...
        'selectedFile', '', 'screenFile', '', 'enzyme_spacing', '', ...
        'outputFile', 'figures/first_result/fig_example_response_curves.png', ...
        'nCurves', 6, 'n_E_vals', 120);
    if mod(nargin, 2) ~= 0
        error('plot_response_curves:InvalidArguments', 'Use name-value arguments.');
    end
    for i = 1:2:nargin
        name = lower(strrep(char(varargin{i}), '-', '_'));
        value = varargin{i+1};
        switch name
            case 'mode', cfg.mode = char(value);
            case 'enzyme_spacing', cfg.enzyme_spacing = char(value);
            case {'selectedfile','selected_file'}, cfg.selectedFile = char(value);
            case {'screenfile','screen_file'}, cfg.screenFile = char(value);
            case {'outputfile','output_file'}, cfg.outputFile = char(value);
            case 'ncurves', cfg.nCurves = value;
            case 'n_e_vals', cfg.n_E_vals = value;
            otherwise
                error('plot_response_curves:UnknownOption', 'Unknown option %s', char(varargin{i}));
        end
    end
end
