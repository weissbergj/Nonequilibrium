function Results = run_deterministic_screen(varargin)
%RUN_DETERMINISTIC_SCREEN Stage 1: deterministic parameter screen (no SSA/KL).
%
% Example (legacy, default):
%   run_deterministic_screen("N_runs",500,"n_E_vals",100,"useParallel",true)
%
% Example (paper-inspired dimensionless):
%   run_deterministic_screen("mode","paper_inspired_dimensionless", ...
%       "N_runs",500,"n_E_vals",100,"useParallel",true)
%
% In dimensionless mode, energy columns are named phi_kBT_E0/E1 and
% delta_mu_kBT_E0/E1 and are reported in kBT-like units (RT=1, ATP=1).

cfg = parse_config(varargin{:});
if ~isempty(cfg.randomSeed)
    rng(cfg.randomSeed);
end

mc = model_constants(cfg.mode);
N = cfg.N_runs;
params = sample_parameter_sets(N, cfg.mode);
out = repmat(empty_result(), N, 1);

fprintf('Deterministic screen [%s]: N_runs=%d, n_E_vals=%d, useParallel=%d\n', ...
    mc.mode, N, cfg.n_E_vals, cfg.useParallel);
if strcmp(mc.input_signal_definition, 'total_enzyme_concentration')
    warning('run_deterministic_screen:InputIsTotalEnzyme', ...
        ['Input signal = TOTAL ENZYME concentration (E0/E1 around E*). This is an ', ...
         'analogy and may NOT reproduce the paper''s input concentration c. ', ...
         'See docs/MODEL_AUDIT.md (TODO: substrate_or_alpha preset).']);
end

if cfg.useParallel
    parfor run = 1:N
        out(run) = screen_one_run(run, params(run,:), cfg.n_E_vals, mc, cfg.enzyme_spacing);
    end
else
    for run = 1:N
        out(run) = screen_one_run(run, params(run,:), cfg.n_E_vals, mc, cfg.enzyme_spacing);
    end
end

Results = outputs_to_table(out, params, cfg, mc);
out_dir = fileparts(cfg.outputFile);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
writetable(Results, cfg.outputFile);
fprintf('Wrote deterministic screen to %s (%d rows, %d valid_for_screen)\n', ...
    cfg.outputFile, height(Results), sum(Results.valid_for_screen));
end

function result = screen_one_run(run, k, n_E_vals, mc, enzyme_spacing)
    result = empty_result();
    result.run_id = run;

    try
        E_vals = build_enzyme_sweep(mc.E_min, mc.E_max, n_E_vals, enzyme_spacing);
        ratio = nan(size(E_vals));
        ode_tol_warn = false;
        ode_failed_any = false;

        for i = 1:numel(E_vals)
            [ratio(i), ode_status] = compute_output_ratio(E_vals(i), mc.alpha, mc.beta, mc.ATP, ...
                k(1), k(2), k(3), k(4), k(5), k(6), mc.t_end, mc.epsilon);
            ode_tol_warn = ode_tol_warn || ode_status.ode_tolerance_warning;
            ode_failed_any = ode_failed_any || ode_status.ode_failed;
        end

        star = find_E_star_E0_E1(E_vals, ratio, mc.E0_frac, mc.E1_frac);
        if star.failed
            result.status = 'failed';
            result.failed_E_star = true;
            result.warning_flags = 'failed_E_star:no_finite_ratio';
            return
        end

        m0 = compute_steady_state_metrics(star.E0, mc.alpha, mc.beta, mc.ATP, ...
            k(1), k(2), k(3), k(4), k(5), k(6), mc.t_end, mc.epsilon, mc.R, mc.T);
        m1 = compute_steady_state_metrics(star.E1, mc.alpha, mc.beta, mc.ATP, ...
            k(1), k(2), k(3), k(4), k(5), k(6), mc.t_end, mc.epsilon, mc.R, mc.T);

        result.E_star = star.E_star;
        result.E_star_index = star.idx_star;
        result.E0 = star.E0;
        result.E1 = star.E1;
        result.ratio_E0 = m0.ratio;
        result.ratio_E1 = m1.ratio;
        result.phi_E0 = m0.phi;
        result.phi_E1 = m1.phi;
        result.delta_mu_E0 = m0.dmu;
        result.delta_mu_E1 = m1.dmu;
        result.Jnet_E0 = m0.Jnet;
        result.Jnet_E1 = m1.Jnet;

        result.ode_tolerance_warning = ode_tol_warn || m0.ode_tolerance_warning || m1.ode_tolerance_warning;
        result.ode_failed = ode_failed_any || m0.ode_failed || m1.ode_failed;
        result.failed_E_star = false;
        result.E_star_at_sweep_boundary = star.E_star_at_sweep_boundary;
        result.complex_phi = ~isreal(result.phi_E0) || ~isreal(result.phi_E1);
        result.negative_phi = real(result.phi_E0) < 0 || real(result.phi_E1) < 0;
        result.nan_or_inf = any(~isfinite([result.ratio_E0, result.ratio_E1, ...
            real(result.phi_E0), real(result.phi_E1), result.delta_mu_E0, result.delta_mu_E1, ...
            result.Jnet_E0, result.Jnet_E1]));

        parts = strings(0,1);
        if result.ode_tolerance_warning, parts(end+1) = "ode_tolerance_warning"; end
        if result.ode_failed, parts(end+1) = "ode_failed"; end
        if result.E_star_at_sweep_boundary, parts(end+1) = "E_star_at_sweep_boundary"; end
        if star.small_ratio_dynamic_range, parts(end+1) = "small_ratio_dynamic_range"; end
        if result.complex_phi, parts(end+1) = "complex_phi"; end
        if result.negative_phi, parts(end+1) = "negative_phi"; end
        if result.nan_or_inf, parts(end+1) = "nan_or_inf"; end

        result.valid_for_screen = ~result.failed_E_star && ~result.ode_failed && ...
            ~result.E_star_at_sweep_boundary && ~result.complex_phi && ...
            ~result.negative_phi && ~result.nan_or_inf;
        result.status = 'ok';
        result.warning_flags = strjoin(parts, ';');
    catch ME
        result.status = 'failed';
        result.failed_E_star = true;
        result.warning_flags = ['error:' ME.identifier ':' ME.message];
    end
end

function result = empty_result()
    result = struct('run_id', nan, 'E_star', nan, 'E_star_index', nan, 'E0', nan, 'E1', nan, ...
        'ratio_E0', nan, 'ratio_E1', nan, 'phi_E0', nan, 'phi_E1', nan, ...
        'delta_mu_E0', nan, 'delta_mu_E1', nan, 'Jnet_E0', nan, 'Jnet_E1', nan, ...
        'ode_tolerance_warning', false, 'ode_failed', false, 'failed_E_star', true, ...
        'E_star_at_sweep_boundary', false, 'complex_phi', false, 'negative_phi', false, ...
        'nan_or_inf', true, 'valid_for_screen', false, 'status', 'not_run', 'warning_flags', '');
end

function T = outputs_to_table(out, params, cfg, mc)
    N = numel(out);
    is_dimensionless = strcmp(mc.mode, 'paper_inspired_dimensionless');

    T = table();
    T.run_id = [out.run_id]';
    T.k1 = params(:,1);
    T.k_minus1 = params(:,2);
    T.k2 = params(:,3);
    T.k_minus2 = params(:,4);
    T.k3 = params(:,5);
    T.k_minus3 = params(:,6);
    T.E_star = [out.E_star]';
    T.E_star_index = [out.E_star_index]';
    T.E0 = [out.E0]';
    T.E1 = [out.E1]';
    T.ratio_E0 = [out.ratio_E0]';
    T.ratio_E1 = [out.ratio_E1]';

    phi_E0 = [out.phi_E0]';
    phi_E1 = [out.phi_E1]';
    dmu_E0 = [out.delta_mu_E0]';
    dmu_E1 = [out.delta_mu_E1]';

    if is_dimensionless
        T.phi_kBT_E0 = phi_E0;
        T.phi_kBT_E1 = phi_E1;
        T.delta_mu_kBT_E0 = dmu_E0;
        T.delta_mu_kBT_E1 = dmu_E1;
    else
        T.phi_E0 = phi_E0;
        T.phi_E1 = phi_E1;
        T.delta_mu_E0 = dmu_E0;
        T.delta_mu_E1 = dmu_E1;
    end

    T.Jnet_E0 = [out.Jnet_E0]';
    T.Jnet_E1 = [out.Jnet_E1]';
    % Canonical aggregate columns used by selection and plotting in both modes.
    T.average_phi = 0.5*(phi_E0 + phi_E1);
    T.average_delta_mu = 0.5*(dmu_E0 + dmu_E1);

    T.ode_tolerance_warning = [out.ode_tolerance_warning]';
    T.ode_failed = [out.ode_failed]';
    T.failed_E_star = [out.failed_E_star]';
    T.E_star_at_sweep_boundary = [out.E_star_at_sweep_boundary]';
    T.complex_phi = [out.complex_phi]';
    T.negative_phi = [out.negative_phi]';
    T.nan_or_inf = [out.nan_or_inf]';
    T.valid_for_screen = [out.valid_for_screen]';
    T.status = string({out.status}');
    T.warning_flags = string({out.warning_flags}');

    T.mode_setting = repmat(string(mc.mode), N, 1);
    T.energy_units = repmat(string(mc.energy_units), N, 1);
    T.enzyme_spacing_setting = repmat(string(cfg.enzyme_spacing), N, 1);
    T.N_runs_setting = repmat(cfg.N_runs, N, 1);
    T.n_E_vals = repmat(cfg.n_E_vals, N, 1);
    T.n_E_vals_setting = repmat(cfg.n_E_vals, N, 1);
    T.useParallel_setting = repmat(cfg.useParallel, N, 1);
    T.randomSeed_setting = repmat(numeric_or_nan(cfg.randomSeed), N, 1);
    T.alpha = repmat(mc.alpha, N, 1);
    T.beta = repmat(mc.beta, N, 1);
    T.ATP = repmat(mc.ATP, N, 1);
    T.RT = repmat(mc.R*mc.T, N, 1);
    T.t_end = repmat(mc.t_end, N, 1);
    T.epsilon = repmat(mc.epsilon, N, 1);
    T.E_min = repmat(mc.E_min, N, 1);
    T.E_max = repmat(mc.E_max, N, 1);
    T.E0_frac = repmat(mc.E0_frac, N, 1);
    T.E1_frac = repmat(mc.E1_frac, N, 1);
    T.input_signal_definition = repmat(string(mc.input_signal_definition), N, 1);
end

function cfg = parse_config(varargin)
    cfg = struct('mode', 'legacy', 'N_runs', 500, 'n_E_vals', 100, 'useParallel', true, ...
        'enzyme_spacing', 'linear', ...
        'outputFile', 'results/first_result/results_deterministic_screen.csv', 'randomSeed', []);
    if mod(nargin, 2) ~= 0
        error('run_deterministic_screen:InvalidArguments', 'Use name-value arguments.');
    end
    for i = 1:2:nargin
        name = normalize_name(varargin{i});
        value = varargin{i+1};
        switch name
            case 'mode', cfg.mode = char(value);
            case 'n_runs', cfg.N_runs = value;
            case 'n_e_vals', cfg.n_E_vals = value;
            case 'enzyme_spacing', cfg.enzyme_spacing = lower(char(value));
            case 'useparallel', cfg.useParallel = logical(value);
            case 'outputfile', cfg.outputFile = char(value);
            case 'randomseed', cfg.randomSeed = value;
            otherwise
                error('run_deterministic_screen:UnknownOption', 'Unknown option %s', char(varargin{i}));
        end
    end
    cfg.N_runs = validate_posint(cfg.N_runs, 'N_runs');
    cfg.n_E_vals = validate_posint(cfg.n_E_vals, 'n_E_vals');
    if ~ismember(cfg.enzyme_spacing, {'linear','log'})
        error('run_deterministic_screen:UnknownSpacing', ...
            'enzyme_spacing must be "linear" or "log" (got "%s").', cfg.enzyme_spacing);
    end
end

function name = normalize_name(value)
    name = lower(strrep(char(value), '-', '_'));
end

function value = validate_posint(value, name)
    validateattributes(value, {'numeric'}, {'scalar','integer','positive'}, ...
        'run_deterministic_screen', name);
    value = double(value);
end

function value = numeric_or_nan(value)
    if isempty(value), value = nan; end
end
