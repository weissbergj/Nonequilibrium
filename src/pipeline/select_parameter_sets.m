function Selected = select_parameter_sets(varargin)
%SELECT_PARAMETER_SETS Stage 2: stratified selection for SSA/KL follow-up.

cfg = parse_config(varargin{:});
if ~isempty(cfg.randomSeed)
    rng(cfg.randomSeed);
end

T = readtable(cfg.inputFile);
valid = valid_screen_rows(T);
V = T(valid,:);

if isempty(V)
    error('select_parameter_sets:NoValidRows', 'No valid deterministic-screen rows found.');
end

switch cfg.stratifyBy
    case 'average_phi'
        metric = real(V.average_phi);
        labels = ["low_energy"; "medium_energy"; "high_energy"];
    case 'average_delta_mu'
        metric = V.average_delta_mu;
        labels = ["low_delta_mu"; "medium_delta_mu"; "high_delta_mu"];
    otherwise
        error('select_parameter_sets:UnknownStratifier', 'Unsupported stratifyBy value.');
end

target_total = min(ceil(cfg.N_select * cfg.oversampleFactor), height(V));
bin_id = stratify_tertiles(metric);
target_per_bin = floor(target_total/3) * ones(3,1);
target_per_bin(1:rem(target_total,3)) = target_per_bin(1:rem(target_total,3)) + 1;

selected_indices = [];
selection_labels = strings(0,1);
for b = 1:3
    candidates = find(bin_id == b);
    if isempty(candidates), continue; end
    take = min(target_per_bin(b), numel(candidates));
    order = candidates(randperm(numel(candidates), take));
    selected_indices = [selected_indices; order(:)]; %#ok<AGROW>
    selection_labels = [selection_labels; repmat(labels(b), take, 1)]; %#ok<AGROW>
end

if numel(selected_indices) < target_total
    remaining = setdiff((1:height(V))', selected_indices);
    take = min(target_total - numel(selected_indices), numel(remaining));
    if take > 0
        extra = remaining(randperm(numel(remaining), take));
        selected_indices = [selected_indices; extra(:)];
        selection_labels = [selection_labels; labels(bin_id(extra))];
    end
end

Selected = V(selected_indices,:);
Selected.selection_bin = selection_labels;
Selected.selection_metric = metric(selected_indices);
Selected.selection_stratifyBy = repmat(string(cfg.stratifyBy), height(Selected), 1);
Selected.selection_source_file = repmat(string(cfg.inputFile), height(Selected), 1);
Selected.selection_randomSeed = repmat(numeric_or_nan(cfg.randomSeed), height(Selected), 1);
Selected.selection_target_total = repmat(cfg.N_select, height(Selected), 1);
Selected.selection_oversample_factor = repmat(cfg.oversampleFactor, height(Selected), 1);

out_dir = fileparts(cfg.outputFile);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
writetable(Selected, cfg.outputFile);
fprintf('Selected %d rows (target=%d, oversample=%g) from %d valid rows. Wrote %s\n', ...
    height(Selected), cfg.N_select, cfg.oversampleFactor, height(V), cfg.outputFile);
end

function valid = valid_screen_rows(T)
    if ismember('valid_for_screen', T.Properties.VariableNames)
        valid = logical(T.valid_for_screen);
        return
    end

    valid = true(height(T), 1);
    if ismember('failed_E_star', T.Properties.VariableNames)
        valid = valid & ~logical(T.failed_E_star);
    end
    if ismember('status', T.Properties.VariableNames)
        valid = valid & strcmp(string(T.status), "ok");
    end
    if ismember('ode_failed', T.Properties.VariableNames)
        valid = valid & ~logical(T.ode_failed);
    end
    if ismember('E_star_at_sweep_boundary', T.Properties.VariableNames)
        valid = valid & ~logical(T.E_star_at_sweep_boundary);
    end
    if ismember('complex_phi', T.Properties.VariableNames)
        valid = valid & ~logical(T.complex_phi);
    end
    if ismember('negative_phi', T.Properties.VariableNames)
        valid = valid & ~logical(T.negative_phi);
    end
    if ismember('nan_or_inf', T.Properties.VariableNames)
        valid = valid & ~logical(T.nan_or_inf);
    end

    required = {'E_star','E0','E1','phi_E0','phi_E1','delta_mu_E0','delta_mu_E1','Jnet_E0','Jnet_E1'};
    for i = 1:numel(required)
        if ismember(required{i}, T.Properties.VariableNames)
            vals = T.(required{i});
            if ~isreal(vals)
                valid = valid & (imag(vals) == 0);
                vals = real(vals);
            end
            valid = valid & isfinite(vals);
        else
            valid = false(height(T), 1);
        end
    end

    if ismember('phi_E0', T.Properties.VariableNames)
        valid = valid & real(T.phi_E0) >= 0 & real(T.phi_E1) >= 0;
    end
end

function bin_id = stratify_tertiles(metric)
    [~, order] = sort(metric);
    bin_id = zeros(size(metric));
    n = numel(metric);
    cut1 = floor(n/3);
    cut2 = floor(2*n/3);
    bin_id(order(1:max(cut1,1))) = 1;
    bin_id(order(max(cut1,1)+1:max(cut2,cut1+1))) = 2;
    bin_id(order(max(cut2,cut1+1)+1:end)) = 3;
    bin_id(bin_id == 0) = 3;
end

function cfg = parse_config(varargin)
    cfg = struct('inputFile', 'results/first_result/results_deterministic_screen.csv', ...
        'outputFile', 'results/first_result/selected_parameter_sets.csv', ...
        'N_select', 50, 'oversampleFactor', 1.25, ...
        'stratifyBy', 'average_phi', 'randomSeed', []);
    if mod(nargin, 2) ~= 0
        error('select_parameter_sets:InvalidArguments', 'Use name-value arguments.');
    end
    for i = 1:2:nargin
        name = normalize_name(varargin{i});
        value = varargin{i+1};
        switch name
            case 'inputfile', cfg.inputFile = char(value);
            case 'outputfile', cfg.outputFile = char(value);
            case 'n_select', cfg.N_select = value;
            case 'oversamplefactor', cfg.oversampleFactor = value;
            case 'stratifyby', cfg.stratifyBy = lower(char(value));
            case 'randomseed', cfg.randomSeed = value;
            otherwise
                error('select_parameter_sets:UnknownOption', 'Unknown option %s', char(varargin{i}));
        end
    end
    cfg.N_select = validate_posint(cfg.N_select, 'N_select');
    validateattributes(cfg.oversampleFactor, {'numeric'}, {'scalar','positive','finite'}, ...
        'select_parameter_sets', 'oversampleFactor');
end

function name = normalize_name(value)
    name = lower(strrep(char(value), '-', '_'));
    switch name
        case {'n_select','nselected','total','n_total'}, name = 'n_select';
        case {'oversamplefactor','oversample_factor'}, name = 'oversamplefactor';
        case {'stratifyby','stratify_by','metric'}, name = 'stratifyby';
    end
end

function value = validate_posint(value, name)
    validateattributes(value, {'numeric'}, {'scalar','integer','positive'}, ...
        'select_parameter_sets', name);
    value = double(value);
end

function value = numeric_or_nan(value)
    if isempty(value), value = nan; end
end
