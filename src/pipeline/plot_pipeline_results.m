function plot_pipeline_results(varargin)
%PLOT_PIPELINE_RESULTS Stage 4: standard KL-vs-energy figures.

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

note = ['Input signal: total enzyme concentration (E0/E1). ', ...
    'Primary points: valid_for_kl only. See docs/RESULT_INTERPRETATION.md.'];

if ~isempty(S) && all(ismember({'average_phi','KL_E0_E1'}, S.Properties.VariableNames))
    make_valid_scatter(S, S_valid, valid_mask, 'average_phi', 'KL_E0_E1', ...
        'Average phi (energy dissipation)', 'KL(E0 || E1)', ...
        'KL vs average phi', note, fullfile(cfg.outputDir, 'fig_KL_vs_phi.png'));
end

if ~isempty(S) && all(ismember({'average_delta_mu','KL_E0_E1'}, S.Properties.VariableNames))
    make_valid_scatter(S, S_valid, valid_mask, 'average_delta_mu', 'KL_E0_E1', ...
        'Average delta mu (J/mol)', 'KL(E0 || E1)', ...
        'KL vs average delta mu', note, fullfile(cfg.outputDir, 'fig_KL_vs_delta_mu.png'));
end

if ~isempty(S) && all(ismember({'average_phi','symmetric_KL'}, S.Properties.VariableNames))
    make_valid_scatter(S, S_valid, valid_mask, 'average_phi', 'symmetric_KL', ...
        'Average phi (energy dissipation)', 'Symmetric KL', ...
        'Symmetric KL vs average phi', note, fullfile(cfg.outputDir, 'fig_symmetric_KL_vs_average_phi.png'));
end

energy_source = choose_plot_source(S, D, 'average_phi');
if ~isempty(energy_source)
    make_histogram(real(energy_source), 'Average phi', ...
        'Deterministic screen: average phi distribution', fullfile(cfg.outputDir, 'fig_energy_distribution.png'));
end

delta_mu_source = choose_plot_source(S, D, 'average_delta_mu');
if ~isempty(delta_mu_source)
    make_histogram(delta_mu_source, 'Average delta mu (J/mol)', ...
        'Deterministic screen: average delta mu distribution', fullfile(cfg.outputDir, 'fig_delta_mu_distribution.png'));
end

if exist(cfg.exampleDistributionFile, 'file')
    E = readtable(cfg.exampleDistributionFile);
    if all(ismember({'condition','output_ratio'}, E.Properties.VariableNames))
        make_example_histograms(E, fullfile(cfg.outputDir, 'fig_example_E0_E1_histograms.png'));
    end
else
    fprintf('No example distribution file; skipping fig_example_E0_E1_histograms.png\n');
end

if ~isempty(S)
    make_validity_summary(S, valid_mask, fullfile(cfg.outputDir, 'fig_validity_summary.png'));
end

fprintf('Figures written to %s\n', cfg.outputDir);
end

function [S_valid, valid_mask] = filter_valid_kl_rows(S)
    if isempty(S), S_valid = S; valid_mask = []; return; end
    if ismember('valid_for_kl', S.Properties.VariableNames)
        valid_mask = logical(S.valid_for_kl);
    else
        valid_mask = isfinite(S.symmetric_KL);
    end
    S_valid = S(valid_mask,:);
    fprintf('Using %d/%d valid_for_kl rows for primary KL trend plots.\n', height(S_valid), height(S));
end

function make_valid_scatter(S_all, S_valid, valid_mask, xcol, ycol, xlab, ylab, title_str, note, output_file)
    fig = figure('Visible','off','Position',[100 100 720 520]);
    if ~isempty(S_all)
        bg = true(height(S_all),1);
        if ~isempty(valid_mask), bg = ~valid_mask; end
        bx = S_all.(xcol); by = S_all.(ycol);
        if ~isreal(bx), bx = real(bx); end
        ok = bg & isfinite(bx) & isfinite(by);
        if any(ok)
            scatter(bx(ok), by(ok), 16, [0.75 0.75 0.75], 'filled', 'MarkerFaceAlpha', 0.25);
            hold on;
        end
    end
    if ~isempty(S_valid)
        x = S_valid.(xcol); y = S_valid.(ycol);
        if ~isreal(x), x = real(x); end
        ok = isfinite(x) & isfinite(y);
        if any(ok)
            scatter(x(ok), y(ok), 28, [0.12 0.45 0.75], 'filled', 'MarkerFaceAlpha', 0.85);
        end
    end
    xlabel(xlab); ylabel(ylab); title(title_str);
    subtitle(note, 'FontSize', 8);
    grid on; box on;
    saveas(fig, output_file); close(fig);
end

function make_histogram(values, xlab, title_str, output_file)
    values = values(isfinite(values));
    fig = figure('Visible','off','Position',[100 100 640 480]);
    histogram(values, min(40, max(10, round(sqrt(numel(values))))));
    xlabel(xlab); ylabel('Count'); title(title_str);
    grid on; box on;
    saveas(fig, output_file); close(fig);
end

function make_example_histograms(E, output_file)
    condition = string(E.condition);
    fig = figure('Visible','off','Position',[100 100 720 520]);
    histogram(E.output_ratio(condition == "E0"), 25, 'Normalization','pdf', 'FaceAlpha',0.55);
    hold on;
    histogram(E.output_ratio(condition == "E1"), 25, 'Normalization','pdf', 'FaceAlpha',0.55);
    xlabel('Stochastic output P/(P+S)'); ylabel('Probability density');
    title('Example E0 vs E1 output distributions (first selected row)');
    legend('E0','E1'); grid on; box on;
    saveas(fig, output_file); close(fig);
end

function make_validity_summary(S, valid_mask, output_file)
    labels = {'valid_for_kl','max_steps_hit','sparse_bins','degenerate_hist','insufficient_bins'};
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
    fig = figure('Visible','off','Position',[100 100 640 480]);
    bar(categorical(labels), counts);
    ylabel('Row count'); title('SSA/KL validity summary');
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
    cfg = struct('screenFile', 'results/first_result/results_deterministic_screen.csv', ...
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
            case {'screenfile','screen_file'}, cfg.screenFile = char(value);
            case {'ssafile','ssa_file'}, cfg.ssaFile = char(value);
            case {'exampledistributionfile','example_distribution_file'}, cfg.exampleDistributionFile = char(value);
            case {'outputdir','output_dir'}, cfg.outputDir = char(value);
            otherwise
                error('plot_pipeline_results:UnknownOption', 'Unknown option %s', char(varargin{i}));
        end
    end
end
