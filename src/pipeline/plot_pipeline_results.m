function plot_pipeline_results(varargin)
%PLOT_PIPELINE_RESULTS Stage 4: KL/JS-vs-energy figures.
%
% Works for both modes. In paper_inspired_dimensionless mode the energy axes
% are labelled in kBT units and the headline plot is JS divergence (bounded,
% does not suffer the raw-KL epsilon saturation).
%
% Optional "mode" argument forces the axis labelling; otherwise it is inferred
% from the energy_units / mode_setting columns in the SSA results file.

cfg = parse_config(varargin{:});
if ~exist(cfg.outputDir, 'dir')
    mkdir(cfg.outputDir);
end

if exist(cfg.ssaFile, 'file')
    S = readtable(cfg.ssaFile);
    [S_valid, valid_mask] = filter_valid_kl_rows(S);
else
    S = table(); S_valid = table(); valid_mask = [];
    warning('plot_pipeline_results:MissingSSAFile', 'SSA/KL file not found: %s', cfg.ssaFile);
end

if exist(cfg.screenFile, 'file')
    D = readtable(cfg.screenFile);
else
    D = table();
    warning('plot_pipeline_results:MissingScreenFile', 'Screen file not found: %s', cfg.screenFile);
end

mode = resolve_mode(cfg.mode, S, D);
is_dimensionless = strcmp(mode, 'paper_inspired_dimensionless');
if is_dimensionless
    phi_lab = 'Average phi (kBT)';
    dmu_lab = 'Average delta mu (kBT)';
    note = 'Enzyme-model analogy. Blue = valid for KL. See docs/MODEL_AUDIT.md.';
else
    phi_lab = 'Average phi (energy dissipation)';
    dmu_lab = 'Average delta mu (J/mol)';
    note = 'Input = total enzyme. Blue = valid for KL. See docs/RESULT_INTERPRETATION.md.';
end

% ---- Headline metric: JS divergence (bounded) ----
if has_cols(S, {'average_phi','JS_divergence'})
    make_valid_scatter(S, S_valid, valid_mask, 'average_phi', 'JS_divergence', ...
        phi_lab, 'Jensen-Shannon divergence (nats, 0..ln2)', ...
        'JS divergence vs energy dissipation: corrected enzyme model', note, ...
        fullfile(cfg.outputDir, 'fig_JS_vs_phi.png'));
end
if has_cols(S, {'average_delta_mu','JS_divergence'})
    make_valid_scatter(S, S_valid, valid_mask, 'average_delta_mu', 'JS_divergence', ...
        dmu_lab, 'Jensen-Shannon divergence (nats, 0..ln2)', ...
        'JS divergence vs average delta mu', note, fullfile(cfg.outputDir, 'fig_JS_vs_delta_mu.png'));
end

% ---- log1p(symmetric KL): compressed but unbounded ----
if has_cols(S, {'average_phi','log1p_symmetric_KL'})
    make_valid_scatter(S, S_valid, valid_mask, 'average_phi', 'log1p_symmetric_KL', ...
        phi_lab, 'log(1 + symmetric KL)', ...
        'log1p symmetric KL vs average phi', note, fullfile(cfg.outputDir, 'fig_log1p_symmetric_KL_vs_phi.png'));
end

% ---- Raw symmetric KL: kept for comparison (may saturate near epsilon ceiling) ----
if has_cols(S, {'average_phi','symmetric_KL'})
    make_valid_scatter(S, S_valid, valid_mask, 'average_phi', 'symmetric_KL', ...
        phi_lab, 'Symmetric KL (raw; can saturate)', ...
        'Raw symmetric KL vs average phi (comparison)', note, fullfile(cfg.outputDir, 'fig_symmetric_KL_vs_average_phi.png'));
end

% ---- Raw directional KL (legacy continuity) ----
if has_cols(S, {'average_phi','KL_E0_E1'})
    make_valid_scatter(S, S_valid, valid_mask, 'average_phi', 'KL_E0_E1', ...
        phi_lab, 'KL(E0 || E1) (raw; can saturate)', ...
        'Raw KL vs average phi', note, fullfile(cfg.outputDir, 'fig_KL_vs_phi.png'));
end
if has_cols(S, {'average_delta_mu','KL_E0_E1'})
    make_valid_scatter(S, S_valid, valid_mask, 'average_delta_mu', 'KL_E0_E1', ...
        dmu_lab, 'KL(E0 || E1) (raw; can saturate)', ...
        'Raw KL vs average delta mu', note, fullfile(cfg.outputDir, 'fig_KL_vs_delta_mu.png'));
end

% ---- Screen distributions ----
energy_source = choose_plot_source(S, D, 'average_phi');
if ~isempty(energy_source)
    make_histogram(real(energy_source), phi_lab, ...
        'Deterministic screen: average phi distribution', fullfile(cfg.outputDir, 'fig_energy_distribution.png'));
end
delta_mu_source = choose_plot_source(S, D, 'average_delta_mu');
if ~isempty(delta_mu_source)
    make_histogram(delta_mu_source, dmu_lab, ...
        'Deterministic screen: average delta mu distribution', fullfile(cfg.outputDir, 'fig_delta_mu_distribution.png'));
end

% ---- Example output distributions ----
if exist(cfg.exampleDistributionFile, 'file')
    E = readtable(cfg.exampleDistributionFile);
    if has_cols(E, {'condition','output_ratio'})
        make_example_histograms(E, fullfile(cfg.outputDir, 'fig_example_E0_E1_histograms.png'));
        if ismember('energy_level', E.Properties.VariableNames)
            make_example_by_energy(E, fullfile(cfg.outputDir, 'fig_example_histograms_by_energy.png'));
        end
    end
else
    fprintf('No example distribution file; skipping example histograms\n');
end

if ~isempty(S)
    make_validity_summary(S, valid_mask, fullfile(cfg.outputDir, 'fig_validity_summary.png'));
end

fprintf('Figures written to %s\n', cfg.outputDir);
end

function tf = has_cols(T, cols)
    tf = ~isempty(T) && all(ismember(cols, T.Properties.VariableNames));
end

function mode = resolve_mode(forced, S, D)
    if ~isempty(forced)
        mode = char(forced);
        if strcmpi(mode, 'dimensionless'), mode = 'paper_inspired_dimensionless'; end
        return
    end
    mode = 'legacy';
    for T = {S, D}
        Tbl = T{1};
        if isempty(Tbl), continue; end
        if ismember('mode_setting', Tbl.Properties.VariableNames) && height(Tbl) > 0
            mode = char(string(Tbl.mode_setting(1)));
            return
        end
        if ismember('energy_units', Tbl.Properties.VariableNames) && height(Tbl) > 0
            if strcmpi(string(Tbl.energy_units(1)), 'kBT')
                mode = 'paper_inspired_dimensionless';
                return
            end
        end
    end
end

function [S_valid, valid_mask] = filter_valid_kl_rows(S)
    if isempty(S), S_valid = S; valid_mask = []; return; end
    if ismember('valid_for_kl', S.Properties.VariableNames)
        valid_mask = logical(S.valid_for_kl);
    else
        valid_mask = isfinite(S.symmetric_KL);
    end
    S_valid = S(valid_mask,:);
    fprintf('Using %d/%d valid_for_kl rows for primary trend plots.\n', height(S_valid), height(S));
end

function make_valid_scatter(S_all, S_valid, valid_mask, xcol, ycol, xlab, ylab, title_str, note, output_file)
    fig = figure('Visible','off','Position',[100 100 720 520]);
    if ~isempty(S_all)
        bg = true(height(S_all),1);
        if ~isempty(valid_mask), bg = ~valid_mask; end
        bx = S_all.(xcol); by = S_all.(ycol);
        if ~isreal(bx), bx = real(bx); end
        if ~isreal(by), by = real(by); end
        ok = bg & isfinite(bx) & isfinite(by);
        if any(ok)
            scatter(bx(ok), by(ok), 16, [0.75 0.75 0.75], 'filled', 'MarkerFaceAlpha', 0.25);
            hold on;
        end
    end
    if ~isempty(S_valid)
        x = S_valid.(xcol); y = S_valid.(ycol);
        if ~isreal(x), x = real(x); end
        if ~isreal(y), y = real(y); end
        ok = isfinite(x) & isfinite(y);
        if any(ok)
            scatter(x(ok), y(ok), 28, [0.12 0.45 0.75], 'filled', 'MarkerFaceAlpha', 0.85);
        end
    end
    xlabel(xlab); ylabel(ylab);
    title(title_str, 'FontSize', 11, 'Interpreter', 'none');
    subtitle(note, 'FontSize', 8, 'Interpreter', 'none');
    grid on; box on;
    saveas(fig, output_file); close(fig);
end

function make_histogram(values, xlab, title_str, output_file)
    values = values(isfinite(values));
    fig = figure('Visible','off','Position',[100 100 640 480]);
    histogram(values, min(40, max(10, round(sqrt(max(1,numel(values)))))));
    xlabel(xlab); ylabel('Count'); title(title_str);
    grid on; box on;
    saveas(fig, output_file); close(fig);
end

function make_example_histograms(E, output_file)
    if ismember('energy_level', E.Properties.VariableNames)
        levels = unique(string(E.energy_level), 'stable');
        pick = pick_one_level(levels);
        E = E(string(E.energy_level) == pick, :);
        title_suffix = char(" (" + pick + " energy example)");
    else
        title_suffix = ' (first selected row)';
    end
    condition = string(E.condition);
    fig = figure('Visible','off','Position',[100 100 720 520]);
    histogram(E.output_ratio(condition == "E0"), 25, 'Normalization','pdf', 'FaceAlpha',0.55);
    hold on;
    histogram(E.output_ratio(condition == "E1"), 25, 'Normalization','pdf', 'FaceAlpha',0.55);
    xlabel('Stochastic output P/(P+S)'); ylabel('Probability density');
    title(['Example E0 vs E1 output distributions' title_suffix]);
    legend('E0','E1'); grid on; box on;
    saveas(fig, output_file); close(fig);
end

function make_example_by_energy(E, output_file)
    order = ["low","medium","high"];
    present = unique(string(E.energy_level), 'stable');
    levels = order(ismember(order, present));
    extra = present(~ismember(present, order));
    levels = [levels, reshape(extra, 1, [])];
    n = numel(levels);
    if n == 0, return; end

    fig = figure('Visible','off','Position',[100 100 360*n 460]);
    for i = 1:n
        subplot(1, n, i);
        sub = E(string(E.energy_level) == levels(i), :);
        cond = string(sub.condition);
        histogram(sub.output_ratio(cond == "E0"), 20, 'Normalization','pdf', 'FaceAlpha',0.55);
        hold on;
        histogram(sub.output_ratio(cond == "E1"), 20, 'Normalization','pdf', 'FaceAlpha',0.55);
        xlabel('P/(P+S)'); ylabel('pdf');
        title(char(levels(i) + " energy"));
        if i == 1, legend('E0','E1','Location','best'); end
        grid on; box on;
    end
    sgtitle('Example E0 vs E1 output distributions by energy level', 'FontSize', 11);
    saveas(fig, output_file); close(fig);
end

function pick = pick_one_level(levels)
    if any(levels == "medium")
        pick = "medium";
    else
        pick = levels(1);
    end
end

function make_validity_summary(S, valid_mask, output_file)
    labels = {'valid KL','max steps hit','sparse bins','degenerate histogram','insufficient bins'};
    counts = zeros(1, numel(labels));
    counts(1) = sum(valid_mask);
    if ismember('ssa_warning_flags', S.Properties.VariableNames)
        flags = string(S.ssa_warning_flags);
        counts(2) = sum(contains(flags, 'max_steps_hit'));
        counts(3) = sum(contains(flags, 'sparse_bins'));
        counts(4) = sum(contains(flags, 'degenerate_histogram'));
        counts(5) = sum(contains(flags, 'insufficient_samples_for_bins'));
    else
        if ismember('max_steps_hits_E0', S.Properties.VariableNames)
            counts(2) = sum(S.max_steps_hits_E0 > 0 | S.max_steps_hits_E1 > 0);
        end
        if ismember('KL_sparse_bins', S.Properties.VariableNames)
            counts(3) = sum(S.KL_sparse_bins);
        end
        if ismember('KL_degenerate_histogram', S.Properties.VariableNames)
            counts(4) = sum(S.KL_degenerate_histogram);
        end
        if ismember('insufficient_samples_for_bins', S.Properties.VariableNames)
            counts(5) = sum(S.insufficient_samples_for_bins);
        end
    end
    fig = figure('Visible','off','Position',[100 100 720 520]);
    bar(1:numel(labels), counts, 0.6, 'FaceColor', [0.12 0.45 0.75]);
    set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels);
    xtickangle(20);
    ylabel('Row count'); title('SSA / KL validity summary', 'FontSize', 12);
    grid on; box on;
    saveas(fig, output_file); close(fig);
end

function values = choose_plot_source(S, D, column_name)
    values = [];
    if ~isempty(S) && ismember(column_name, S.Properties.VariableNames)
        values = S.(column_name);
    elseif ~isempty(D) && ismember(column_name, D.Properties.VariableNames)
        values = D.(column_name);
    end
end

function cfg = parse_config(varargin)
    cfg = struct('mode', '', ...
        'screenFile', 'results/first_result/results_deterministic_screen.csv', ...
        'ssaFile', 'results/first_result/results_ssa_kl.csv', ...
        'exampleDistributionFile', 'results/first_result/ssa_example_E0_E1_distributions.csv', ...
        'outputDir', 'figures/first_result');
    if mod(nargin, 2) ~= 0
        error('plot_pipeline_results:InvalidArguments', 'Use name-value arguments.');
    end
    for i = 1:2:nargin
        name = lower(strrep(char(varargin{i}), '-', '_'));
        value = varargin{i+1};
        switch name
            case 'mode', cfg.mode = char(value);
            case {'screenfile','screen_file'}, cfg.screenFile = char(value);
            case {'ssafile','ssa_file'}, cfg.ssaFile = char(value);
            case {'exampledistributionfile','example_distribution_file'}, cfg.exampleDistributionFile = char(value);
            case {'outputdir','output_dir'}, cfg.outputDir = char(value);
            otherwise
                error('plot_pipeline_results:UnknownOption', 'Unknown option %s', char(varargin{i}));
        end
    end
end
