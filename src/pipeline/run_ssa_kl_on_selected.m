function Results = run_ssa_kl_on_selected(varargin)
%RUN_SSA_KL_ON_SELECTED Stage 3: Gillespie SSA and KL on selected rows only.
%
% Example (legacy):
%   run_ssa_kl_on_selected("inputFile","selected.csv","N_stoch",50)
%
% Example (paper-inspired dimensionless):
%   run_ssa_kl_on_selected("mode","paper_inspired_dimensionless", ...
%       "inputFile","selected.csv","N_stoch",50,"t_end_stoch",25,"nBins",10)
%
% Adds raw KL, symmetric_KL and bounded JS_divergence. JS is the preferred
% interpretive metric because raw KL saturates near -ln(epsilon) when the two
% output histograms do not overlap.

cfg = parse_config(varargin{:});
mc = model_constants(cfg.mode);
if ~isempty(cfg.Vol_factor)
    mc.Vol_factor = cfg.Vol_factor;   % opt-in override (e.g. low-copy-number test)
end
if ~isempty(cfg.max_steps)
    mc.max_steps = cfg.max_steps;     % opt-in override to bound SSA cost
end
if isempty(cfg.epsilon)
    cfg.epsilon = mc.kl_epsilon;
end
if isempty(cfg.t_end_stoch)
    cfg.t_end_stoch = mc.t_end_stoch_default;
end
validateattributes(cfg.t_end_stoch, {'numeric'}, {'scalar','positive','finite'}, ...
    'run_ssa_kl_on_selected', 't_end_stoch');
validateattributes(cfg.epsilon, {'numeric'}, {'scalar','positive','finite'}, ...
    'run_ssa_kl_on_selected', 'epsilon');

if ~isempty(cfg.randomSeed)
    rng(cfg.randomSeed);
end

Selected = readtable(cfg.inputFile);
validate_selected_table(Selected);

N = height(Selected);
out = repmat(empty_result(), N, 1);
batchSize = min(cfg.batchSize, N);

fprintf('SSA/KL on selected sets [%s]: N=%d, N_stoch=%d, t_end_stoch=%g, nBins=%d\n', ...
    mc.mode, N, cfg.N_stoch, cfg.t_end_stoch, cfg.nBins);
warn_if_kl_sampling_is_thin(cfg.N_stoch, cfg.nBins);

for start_idx = 1:batchSize:N
    stop_idx = min(start_idx + batchSize - 1, N);
    batch_indices = start_idx:stop_idx;
    batch_out = repmat(empty_result(), numel(batch_indices), 1);

    if cfg.useParallel
        parfor jj = 1:numel(batch_indices)
            row_idx = batch_indices(jj);
            batch_out(jj) = ssa_one_row(Selected(row_idx,:), row_idx, cfg, mc);
        end
    else
        for jj = 1:numel(batch_indices)
            row_idx = batch_indices(jj);
            batch_out(jj) = ssa_one_row(Selected(row_idx,:), row_idx, cfg, mc);
        end
    end

    out(batch_indices) = batch_out;
end

Results = outputs_to_table(Selected, out, cfg, mc);
out_dir = fileparts(cfg.outputFile);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
writetable(Results, cfg.outputFile);
fprintf('Wrote SSA/KL results to %s (%d rows, %d valid_for_kl)\n', ...
    cfg.outputFile, height(Results), sum(Results.valid_for_kl));

if cfg.saveExampleDistributions && N > 0
    write_example_distributions(out, Selected, cfg);
end
end

function result = ssa_one_row(row, row_idx, cfg, mc)
    result = empty_result();
    result.selected_row = row_idx;
    result.batch_id = ceil(row_idx / cfg.batchSize);
    warning_parts = strings(0,1);

    if ~isempty(cfg.randomSeed)
        rng(cfg.randomSeed + row_idx);
    end

    try
        k1 = row.k1; k_1 = row.k_minus1; k2 = row.k2; k_2 = row.k_minus2;
        k3 = row.k3; k_3 = row.k_minus3;
        E0 = row.E0; E1 = row.E1;

        dist_E0 = nan(cfg.N_stoch, 1);
        dist_E1 = nan(cfg.N_stoch, 1);
        steps_E0 = nan(cfg.N_stoch, 1);
        steps_E1 = nan(cfg.N_stoch, 1);
        hitmax_E0 = false(cfg.N_stoch, 1);
        hitmax_E1 = false(cfg.N_stoch, 1);
        failed_E0 = false(cfg.N_stoch, 1);
        failed_E1 = false(cfg.N_stoch, 1);

        for i = 1:cfg.N_stoch
            [dist_E0(i), steps_E0(i), hitmax_E0(i), failed_E0(i)] = run_gillespie_ssa( ...
                round(E0 * mc.Vol_factor), mc.Vol_factor, mc.alpha, k1, k2, k_3, k_1, k_2, ...
                k3, mc.beta, mc.ATP, cfg.t_end_stoch, mc.max_steps, cfg.epsilon);
            [dist_E1(i), steps_E1(i), hitmax_E1(i), failed_E1(i)] = run_gillespie_ssa( ...
                round(E1 * mc.Vol_factor), mc.Vol_factor, mc.alpha, k1, k2, k_3, k_1, k_2, ...
                k3, mc.beta, mc.ATP, cfg.t_end_stoch, mc.max_steps, cfg.epsilon);
        end

        result.failed_trajectories_E0 = sum(failed_E0);
        result.failed_trajectories_E1 = sum(failed_E1);
        result.max_steps_hits_E0 = sum(hitmax_E0);
        result.max_steps_hits_E1 = sum(hitmax_E1);
        result.mean_steps_E0 = finite_mean(steps_E0);
        result.mean_steps_E1 = finite_mean(steps_E1);

        if any(failed_E0) || any(failed_E1)
            warning_parts(end+1) = "failed_trajectories";
        end
        if any(hitmax_E0) || any(hitmax_E1)
            warning_parts(end+1) = "max_steps_hit";
        end

        kl = compute_kl_divergence(dist_E0, dist_E1, cfg.nBins, cfg.epsilon, cfg.N_stoch);
        warning_parts = [warning_parts, kl.warning_flags]; %#ok<AGROW>

        result.KL_E0_E1 = kl.KL_E0_E1;
        result.KL_E1_E0 = kl.KL_E1_E0;
        result.symmetric_KL = kl.symmetric_KL;
        result.log1p_symmetric_KL = log1p(max(0, kl.symmetric_KL));
        result.JS_divergence = kl.JS_divergence;
        result.KL_degenerate_histogram = kl.KL_degenerate_histogram;
        result.KL_sparse_bins = kl.KL_sparse_bins;
        result.insufficient_samples_for_bins = kl.insufficient_samples_for_bins;
        result.output_mean_E0 = finite_mean(dist_E0);
        result.output_mean_E1 = finite_mean(dist_E1);
        result.output_var_E0 = finite_var(dist_E0);
        result.output_var_E1 = finite_var(dist_E1);
        result.fraction_at_0_E0 = frac_at_boundary(dist_E0, 0, cfg.boundary_tol);
        result.fraction_at_1_E0 = frac_at_boundary(dist_E0, 1, cfg.boundary_tol);
        result.fraction_at_0_E1 = frac_at_boundary(dist_E1, 0, cfg.boundary_tol);
        result.fraction_at_1_E1 = frac_at_boundary(dist_E1, 1, cfg.boundary_tol);
        result.dist_E0 = dist_E0;
        result.dist_E1 = dist_E1;

        finite_real_kl = isreal(kl.KL_E0_E1) && isreal(kl.KL_E1_E0) && ...
            isfinite(kl.KL_E0_E1) && isfinite(kl.KL_E1_E0) && isfinite(kl.symmetric_KL);

        result.valid_for_kl = ~any(failed_E0) && ~any(failed_E1) && ...
            ~any(hitmax_E0) && ~any(hitmax_E1) && ~kl.KL_degenerate_histogram && ...
            ~kl.KL_sparse_bins && ~kl.insufficient_samples_for_bins && finite_real_kl;

        if ~result.valid_for_kl
            warning_parts(end+1) = "invalid_for_kl";
        end

        result.status = 'ok';
        result.warning_flags = strjoin(unique(warning_parts), ';');
    catch ME
        result.status = 'failed';
        result.warning_flags = "error:" + string(ME.identifier) + ":" + string(ME.message);
    end
end

function Results = outputs_to_table(Selected, out, cfg, mc)
    Results = Selected;
    Results.ssa_selected_row = [out.selected_row]';
    Results.ssa_batch_id = [out.batch_id]';
    Results.KL_E0_E1 = [out.KL_E0_E1]';
    Results.KL_E1_E0 = [out.KL_E1_E0]';
    Results.symmetric_KL = [out.symmetric_KL]';
    Results.log1p_symmetric_KL = [out.log1p_symmetric_KL]';
    Results.JS_divergence = [out.JS_divergence]';
    Results.output_mean_E0 = [out.output_mean_E0]';
    Results.output_mean_E1 = [out.output_mean_E1]';
    Results.output_var_E0 = [out.output_var_E0]';
    Results.output_var_E1 = [out.output_var_E1]';
    Results.fraction_at_0_E0 = [out.fraction_at_0_E0]';
    Results.fraction_at_1_E0 = [out.fraction_at_1_E0]';
    Results.fraction_at_0_E1 = [out.fraction_at_0_E1]';
    Results.fraction_at_1_E1 = [out.fraction_at_1_E1]';
    Results.mean_steps_E0 = [out.mean_steps_E0]';
    Results.mean_steps_E1 = [out.mean_steps_E1]';
    Results.failed_trajectories_E0 = [out.failed_trajectories_E0]';
    Results.failed_trajectories_E1 = [out.failed_trajectories_E1]';
    Results.max_steps_hits_E0 = [out.max_steps_hits_E0]';
    Results.max_steps_hits_E1 = [out.max_steps_hits_E1]';
    Results.KL_degenerate_histogram = [out.KL_degenerate_histogram]';
    Results.KL_sparse_bins = [out.KL_sparse_bins]';
    Results.insufficient_samples_for_bins = [out.insufficient_samples_for_bins]';
    Results.valid_for_kl = [out.valid_for_kl]';
    Results.ssa_status = string({out.status}');
    Results.ssa_warning_flags = string({out.warning_flags}');
    Results.mode_setting = repmat(string(mc.mode), height(Results), 1);
    Results.energy_units = repmat(string(mc.energy_units), height(Results), 1);
    Results.N_stoch_setting = repmat(cfg.N_stoch, height(Results), 1);
    Results.t_end_stoch_setting = repmat(cfg.t_end_stoch, height(Results), 1);
    Results.nBins_setting = repmat(cfg.nBins, height(Results), 1);
    Results.epsilon_setting = repmat(cfg.epsilon, height(Results), 1);
    Results.useParallel_ssa_setting = repmat(cfg.useParallel, height(Results), 1);
    Results.batchSize_setting = repmat(cfg.batchSize, height(Results), 1);
    Results.randomSeed_ssa_setting = repmat(numeric_or_nan(cfg.randomSeed), height(Results), 1);
    Results.alpha_ssa = repmat(mc.alpha, height(Results), 1);
    Results.beta_ssa = repmat(mc.beta, height(Results), 1);
    Results.ATP_ssa = repmat(mc.ATP, height(Results), 1);
    Results.Vol_factor = repmat(mc.Vol_factor, height(Results), 1);
    Results.max_steps = repmat(mc.max_steps, height(Results), 1);
end

function write_example_distributions(out, Selected, cfg)
    reps = pick_representative_rows(out, Selected);
    if isempty(reps)
        return
    end

    energy_level = strings(0,1);
    selection_bin = strings(0,1);
    run_id = [];
    selected_row = [];
    condition = strings(0,1);
    sample_id = [];
    output_ratio = [];

    for r = 1:numel(reps)
        idx = reps(r).idx;
        level = reps(r).level;
        bin_label = reps(r).bin_label;
        d0 = out(idx).dist_E0(:);
        d1 = out(idx).dist_E1(:);
        n0 = numel(d0);
        n1 = numel(d1);
        rid = example_run_id(Selected, idx);

        energy_level = [energy_level; repmat(level, n0+n1, 1)]; %#ok<AGROW>
        selection_bin = [selection_bin; repmat(bin_label, n0+n1, 1)]; %#ok<AGROW>
        run_id = [run_id; repmat(rid, n0+n1, 1)]; %#ok<AGROW>
        selected_row = [selected_row; repmat(idx, n0+n1, 1)]; %#ok<AGROW>
        condition = [condition; repmat("E0", n0, 1); repmat("E1", n1, 1)]; %#ok<AGROW>
        sample_id = [sample_id; (1:n0)'; (1:n1)']; %#ok<AGROW>
        output_ratio = [output_ratio; d0; d1]; %#ok<AGROW>
    end

    Example = table(energy_level, selection_bin, run_id, selected_row, condition, sample_id, output_ratio);
    Example.N_stoch_setting = repmat(cfg.N_stoch, height(Example), 1);
    Example.t_end_stoch_setting = repmat(cfg.t_end_stoch, height(Example), 1);

    ex_dir = fileparts(cfg.exampleDistributionFile);
    if ~isempty(ex_dir) && ~exist(ex_dir, 'dir')
        mkdir(ex_dir);
    end
    writetable(Example, cfg.exampleDistributionFile);
    fprintf('Wrote example SSA distributions (%d energy levels) to %s\n', ...
        numel(reps), cfg.exampleDistributionFile);
end

function reps = pick_representative_rows(out, Selected)
    reps = struct('idx', {}, 'level', {}, 'bin_label', {});
    has_bin = ismember('selection_bin', Selected.Properties.VariableNames);
    if has_bin
        bin_strings = string(Selected.selection_bin);
    else
        bin_strings = strings(height(Selected), 1);
    end
    valid = arrayfun(@(o) ~isempty(o.dist_E0) && o.valid_for_kl, out);
    has_dist = arrayfun(@(o) ~isempty(o.dist_E0), out);

    levels = ["low", "medium", "high"];
    used = false(numel(out), 1);
    for L = 1:numel(levels)
        level = levels(L);
        if has_bin
            in_bin = startsWith(bin_strings, level);
        else
            in_bin = true(numel(out), 1);
        end
        idx = find(in_bin(:) & valid(:) & ~used, 1);
        if isempty(idx)
            idx = find(in_bin(:) & has_dist(:) & ~used, 1);
        end
        if isempty(idx) && ~has_bin
            idx = find(has_dist(:) & ~used, 1);
        end
        if ~isempty(idx)
            used(idx) = true;
            if has_bin
                bin_label = bin_strings(idx);
            else
                bin_label = level;
            end
            reps(end+1) = struct('idx', idx, 'level', level, 'bin_label', bin_label); %#ok<AGROW>
        end
    end

    if isempty(reps)
        idx = find(has_dist, 1);
        if ~isempty(idx)
            reps(end+1) = struct('idx', idx, 'level', "first", 'bin_label', "first");
        end
    end
end

function rid = example_run_id(Selected, idx)
    if ismember('run_id', Selected.Properties.VariableNames)
        rid = Selected.run_id(idx);
    else
        rid = idx;
    end
end

function result = empty_result()
    result = struct('selected_row', nan, 'batch_id', nan, ...
        'KL_E0_E1', nan, 'KL_E1_E0', nan, 'symmetric_KL', nan, ...
        'log1p_symmetric_KL', nan, 'JS_divergence', nan, ...
        'output_mean_E0', nan, 'output_mean_E1', nan, ...
        'output_var_E0', nan, 'output_var_E1', nan, ...
        'fraction_at_0_E0', nan, 'fraction_at_1_E0', nan, ...
        'fraction_at_0_E1', nan, 'fraction_at_1_E1', nan, ...
        'mean_steps_E0', nan, 'mean_steps_E1', nan, ...
        'failed_trajectories_E0', nan, 'failed_trajectories_E1', nan, ...
        'max_steps_hits_E0', nan, 'max_steps_hits_E1', nan, ...
        'KL_degenerate_histogram', false, 'KL_sparse_bins', false, ...
        'insufficient_samples_for_bins', false, 'valid_for_kl', false, ...
        'status', 'not_run', 'warning_flags', '', ...
        'dist_E0', [], 'dist_E1', []);
end

function validate_selected_table(T)
    required = {'k1','k_minus1','k2','k_minus2','k3','k_minus3','E0','E1'};
    for i = 1:numel(required)
        if ~ismember(required{i}, T.Properties.VariableNames)
            error('run_ssa_kl_on_selected:MissingColumn', 'Missing column: %s', required{i});
        end
    end
end

function cfg = parse_config(varargin)
    cfg = struct('mode', 'legacy', ...
        'inputFile', 'results/first_result/selected_parameter_sets.csv', ...
        'outputFile', 'results/first_result/results_ssa_kl.csv', ...
        'N_stoch', 50, 't_end_stoch', [], ...
        'nBins', 10, 'epsilon', [], ...
        'useParallel', true, 'batchSize', 25, 'randomSeed', [], ...
        'boundary_tol', 1e-3, 'Vol_factor', [], 'max_steps', [], ...
        'saveExampleDistributions', true, ...
        'exampleDistributionFile', 'results/first_result/ssa_example_E0_E1_distributions.csv');
    if mod(nargin, 2) ~= 0
        error('run_ssa_kl_on_selected:InvalidArguments', 'Use name-value arguments.');
    end
    for i = 1:2:nargin
        name = normalize_name(varargin{i});
        value = varargin{i+1};
        switch name
            case 'mode', cfg.mode = char(value);
            case 'inputfile', cfg.inputFile = char(value);
            case 'outputfile', cfg.outputFile = char(value);
            case 'n_stoch', cfg.N_stoch = value;
            case 't_end_stoch', cfg.t_end_stoch = value;
            case 'nbins', cfg.nBins = value;
            case 'epsilon', cfg.epsilon = value;
            case 'useparallel', cfg.useParallel = logical(value);
            case 'batchsize', cfg.batchSize = value;
            case 'randomseed', cfg.randomSeed = value;
            case 'boundary_tol', cfg.boundary_tol = value;
            case 'vol_factor', cfg.Vol_factor = value;
            case 'max_steps', cfg.max_steps = value;
            case 'saveexampledistributions', cfg.saveExampleDistributions = logical(value);
            case 'exampledistributionfile', cfg.exampleDistributionFile = char(value);
            otherwise
                error('run_ssa_kl_on_selected:UnknownOption', 'Unknown option %s', char(varargin{i}));
        end
    end
    cfg.N_stoch = validate_posint(cfg.N_stoch, 'N_stoch');
    cfg.nBins = validate_posint(cfg.nBins, 'nBins');
    cfg.batchSize = validate_posint(cfg.batchSize, 'batchSize');
end

function warn_if_kl_sampling_is_thin(N_stoch, nBins)
    if nBins > N_stoch/2
        warning('run_ssa_kl_on_selected:BinsExceedHalfSamples', ...
            'nBins=%d > N_stoch/2=%g: histogram KL will be sparse.', nBins, N_stoch/2);
    end
    if N_stoch < 5*nBins
        warning('run_ssa_kl_on_selected:FewSamplesPerBin', ...
            'N_stoch=%d < 5*nBins=%d: KL estimates likely unstable.', N_stoch, 5*nBins);
    end
end

function name = normalize_name(value)
    name = lower(strrep(char(value), '-', '_'));
end

function value = validate_posint(value, name)
    validateattributes(value, {'numeric'}, {'scalar','integer','positive'}, ...
        'run_ssa_kl_on_selected', name);
    value = double(value);
end

function value = finite_mean(x)
    x = x(isfinite(x));
    if isempty(x), value = nan; else, value = mean(x); end
end

function value = finite_var(x)
    x = x(isfinite(x));
    if numel(x) < 2, value = nan; else, value = var(x); end
end

function value = frac_at_boundary(x, boundary, tol)
    x = x(isfinite(x));
    if isempty(x), value = nan; return; end
    if boundary == 0
        value = mean(x <= tol);
    else
        value = mean(x >= 1 - tol);
    end
end

function value = numeric_or_nan(value)
    if isempty(value), value = nan; end
end
