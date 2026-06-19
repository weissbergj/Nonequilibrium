function summary = write_pipeline_summary(screenFile, selectedFile, ssaFile, outputFile, settings)
%WRITE_PIPELINE_SUMMARY Write markdown/text summary of pipeline run counts.

if nargin < 5
    settings = struct();
end

summary = struct();
summary.screen_total = 0;
summary.screen_valid = 0;
summary.selected_total = 0;
summary.ssa_total = 0;
summary.ssa_valid_kl = 0;
summary.ssa_max_steps = 0;
summary.ssa_sparse_bins = 0;
summary.ssa_degenerate = 0;
summary.ssa_insufficient_bins = 0;
summary.ssa_failed_trajectories = 0;
summary.audit = struct();

summary.enzyme_spacing = "";
summary.screen_E_star_boundary = 0;
if exist(screenFile, 'file')
    D = readtable(screenFile);
    summary.screen_total = height(D);
    if ismember('valid_for_screen', D.Properties.VariableNames)
        summary.screen_valid = sum(D.valid_for_screen);
    else
        summary.screen_valid = sum(strcmp(string(D.status), "ok"));
    end
    if ismember('E_star_at_sweep_boundary', D.Properties.VariableNames)
        summary.screen_E_star_boundary = sum(logical(D.E_star_at_sweep_boundary));
    end
    if ismember('enzyme_spacing_setting', D.Properties.VariableNames) && height(D) > 0
        summary.enzyme_spacing = string(D.enzyme_spacing_setting(1));
    end
end

if exist(selectedFile, 'file')
    S = readtable(selectedFile);
    summary.selected_total = height(S);
end

if exist(ssaFile, 'file')
    K = readtable(ssaFile);
    summary.ssa_total = height(K);
    if ismember('valid_for_kl', K.Properties.VariableNames)
        summary.ssa_valid_kl = sum(K.valid_for_kl);
    end
    if ismember('max_steps_hits_E0', K.Properties.VariableNames)
        summary.ssa_max_steps = sum(K.max_steps_hits_E0 > 0 | K.max_steps_hits_E1 > 0);
    end
    if ismember('KL_sparse_bins', K.Properties.VariableNames)
        summary.ssa_sparse_bins = sum(K.KL_sparse_bins);
    end
    if ismember('KL_degenerate_histogram', K.Properties.VariableNames)
        summary.ssa_degenerate = sum(K.KL_degenerate_histogram);
    end
    if ismember('insufficient_samples_for_bins', K.Properties.VariableNames)
        summary.ssa_insufficient_bins = sum(K.insufficient_samples_for_bins);
    end
    if ismember('failed_trajectories_E0', K.Properties.VariableNames)
        summary.ssa_failed_trajectories = sum(K.failed_trajectories_E0 > 0 | K.failed_trajectories_E1 > 0);
    end
    summary.audit = compute_audit(K);
end

lines = strings(0, 1);
lines = append_line(lines, "# Pipeline summary");
lines = append_line(lines, "");
if strlength(summary.enzyme_spacing) > 0
    lines = append_line(lines, sprintf("- Enzyme sweep spacing: **%s**", summary.enzyme_spacing));
end
lines = append_line(lines, sprintf("- Deterministic rows: **%d**", summary.screen_total));
lines = append_line(lines, sprintf("- Valid for screen: **%d**", summary.screen_valid));
lines = append_line(lines, sprintf("- E_star_at_sweep_boundary rows: **%d**", summary.screen_E_star_boundary));
lines = append_line(lines, sprintf("- Selected for SSA: **%d**", summary.selected_total));
lines = append_line(lines, sprintf("- SSA rows: **%d**", summary.ssa_total));
lines = append_line(lines, sprintf("- valid_for_kl: **%d**", summary.ssa_valid_kl));
lines = append_line(lines, sprintf("- max_steps_hit rows: **%d**", summary.ssa_max_steps));
lines = append_line(lines, sprintf("- sparse_bins rows: **%d**", summary.ssa_sparse_bins));
lines = append_line(lines, sprintf("- degenerate_histogram rows: **%d**", summary.ssa_degenerate));
lines = append_line(lines, sprintf("- insufficient_samples_for_bins rows: **%d**", summary.ssa_insufficient_bins));
lines = append_line(lines, sprintf("- failed_trajectories rows: **%d**", summary.ssa_failed_trajectories));
lines = append_line(lines, "");

a = summary.audit;
if isfield(a, 'n_valid') && a.n_valid > 0
    lines = append_line(lines, "## Scientific sanity audit (valid_for_kl rows)");
    lines = append_line(lines, "");

    if isfield(a, 'phi_dmu_max_relerr')
        lines = append_line(lines, "### phi vs delta_mu identity");
        lines = append_line(lines, sprintf("- max relative error of phi vs Jnet*delta_mu: **%.2e** (median %.2e)", ...
            a.phi_dmu_max_relerr, a.phi_dmu_med_relerr));
        lines = append_line(lines, sprintf("- Jnet_E0 min/median/max: %.4g / %.4g / %.4g", a.jnet_min, a.jnet_med, a.jnet_max));
        lines = append_line(lines, sprintf("- fraction of rows with |Jnet-1| < 0.05: **%.2f**", a.jnet_frac_near1));
        lines = append_line(lines, "- Interpretation: phi and delta_mu are NOT duplicated columns; they satisfy the");
        lines = append_line(lines, "  steady-state identity phi = Jnet * delta_mu (entropy production = flux x affinity).");
        lines = append_line(lines, "  Because net cycle flux Jnet ~ alpha = 1 at steady state, phi ~ delta_mu numerically,");
        lines = append_line(lines, "  so JS-vs-phi and JS-vs-delta_mu are nearly the same plot under current settings.");
        lines = append_line(lines, "");
    end

    if isfield(a, 'out_mean_E0')
        lines = append_line(lines, "### Stochastic output distribution");
        lines = append_line(lines, sprintf("- output_mean E0/E1 (median over rows): %.3f / %.3f", a.out_mean_E0, a.out_mean_E1));
        lines = append_line(lines, sprintf("- output_var  E0/E1 (median over rows): %.3f / %.3f", a.out_var_E0, a.out_var_E1));
        lines = append_line(lines, sprintf("- fraction_at_0 E0/E1 (median): %.2f / %.2f", a.frac0_E0, a.frac0_E1));
        lines = append_line(lines, sprintf("- fraction_at_1 E0/E1 (median): %.2f / %.2f", a.frac1_E0, a.frac1_E1));
        lines = append_line(lines, sprintf("- rows where fraction_at_0 + fraction_at_1 > 0.5 (boundary-dominated): **%d/%d**", ...
            a.boundary_dominated, a.n_valid));
        lines = append_line(lines, "- Note: boundary pile-up at 0 and 1 reflects low molecule counts (Vol_factor=1)");
        lines = append_line(lines, "  and/or short t_end_stoch; outputs are near-deterministic per trajectory.");
        lines = append_line(lines, "");
    end

    if isfield(a, 'spearman_phi_js')
        lines = append_line(lines, "### Energy-discrimination correlations (Spearman)");
        lines = append_line(lines, sprintf("- average_phi vs JS_divergence: rho = **%.3f** (n=%d)", a.spearman_phi_js, a.n_valid));
        lines = append_line(lines, sprintf("- average_phi vs log1p_symmetric_KL: rho = **%.3f** (n=%d)", a.spearman_phi_log1pkl, a.n_valid));
        if abs(a.spearman_phi_js) < 0.3
            lines = append_line(lines, "- **The corrected enzyme model does NOT show a strong monotonic");
            lines = append_line(lines, "  energy-discrimination trend under current settings** (|rho| < 0.3).");
        end
        lines = append_line(lines, "");
    end
end

if isfield(settings, 'name') && ~isempty(settings.name)
    lines = append_line(lines, sprintf("Run profile: **%s**", string(settings.name)));
    lines = append_line(lines, "");
end
if isfield(settings, 'notes') && ~isempty(settings.notes)
    note_lines = notes_to_strings(settings.notes);
    for i = 1:numel(note_lines)
        lines = append_line(lines, note_lines(i));
    end
    lines = append_line(lines, "");
end

lines = append_line(lines, "## Caveat");
lines = append_line(lines, "E0/E1 are defined by total enzyme concentration. If the intended input signal is substrate S or influx alpha, change the scientific definition before final interpretation.");
lines = append_line(lines, "");
lines = append_line(lines, "Primary KL trend plots use only rows with valid_for_kl=true.");

out_dir = fileparts(outputFile);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fid = fopen(outputFile, 'w');
if fid < 0
    error('write_pipeline_summary:WriteFailed', 'Could not write %s', outputFile);
end
fprintf(fid, '%s\n', lines);
fclose(fid);

fprintf('\n=== Pipeline summary ===\n');
fprintf('Deterministic rows: %d (valid: %d)\n', summary.screen_total, summary.screen_valid);
fprintf('Selected rows: %d\n', summary.selected_total);
fprintf('SSA rows: %d (valid_for_kl: %d)\n', summary.ssa_total, summary.ssa_valid_kl);
fprintf('Invalid SSA: max_steps=%d, sparse_bins=%d, degenerate=%d\n', ...
    summary.ssa_max_steps, summary.ssa_sparse_bins, summary.ssa_degenerate);
fprintf('Summary written to %s\n', outputFile);

end

function lines = append_line(lines, text)
    lines(end+1, 1) = string(text);
end

function a = compute_audit(K)
    a = struct();
    if isempty(K), a.n_valid = 0; return; end
    if ismember('valid_for_kl', K.Properties.VariableNames)
        v = logical(K.valid_for_kl);
    else
        v = true(height(K), 1);
    end
    Kv = K(v, :);
    a.n_valid = height(Kv);
    if a.n_valid == 0, return; end

    phi0 = get_col(Kv, {'phi_kBT_E0','phi_E0'});
    dmu0 = get_col(Kv, {'delta_mu_kBT_E0','delta_mu_E0'});
    jnet0 = get_col(Kv, {'Jnet_E0'});
    if ~isempty(phi0) && ~isempty(dmu0) && ~isempty(jnet0)
        pred = jnet0 .* dmu0;
        relerr = abs(phi0 - pred) ./ max(abs(phi0), 1e-12);
        relerr = relerr(isfinite(relerr));
        a.phi_dmu_max_relerr = max(relerr);
        a.phi_dmu_med_relerr = median(relerr);
        jf = jnet0(isfinite(jnet0));
        a.jnet_min = min(jf); a.jnet_med = median(jf); a.jnet_max = max(jf);
        a.jnet_frac_near1 = mean(abs(jf - 1) < 0.05);
    end

    om0 = get_col(Kv, {'output_mean_E0'}); om1 = get_col(Kv, {'output_mean_E1'});
    ov0 = get_col(Kv, {'output_var_E0'});  ov1 = get_col(Kv, {'output_var_E1'});
    f00 = get_col(Kv, {'fraction_at_0_E0'}); f10 = get_col(Kv, {'fraction_at_1_E0'});
    f01 = get_col(Kv, {'fraction_at_0_E1'}); f11 = get_col(Kv, {'fraction_at_1_E1'});
    if ~isempty(om0)
        a.out_mean_E0 = nanmedian_local(om0); a.out_mean_E1 = nanmedian_local(om1);
        a.out_var_E0 = nanmedian_local(ov0);  a.out_var_E1 = nanmedian_local(ov1);
    end
    if ~isempty(f00)
        a.frac0_E0 = nanmedian_local(f00); a.frac1_E0 = nanmedian_local(f10);
        a.frac0_E1 = nanmedian_local(f01); a.frac1_E1 = nanmedian_local(f11);
        boundary_mass = f00 + f10;
        a.boundary_dominated = sum(boundary_mass(isfinite(boundary_mass)) > 0.5);
    end

    aphi = get_col(Kv, {'average_phi'});
    js = get_col(Kv, {'JS_divergence'});
    lkl = get_col(Kv, {'log1p_symmetric_KL'});
    if ~isempty(aphi) && ~isempty(js)
        a.spearman_phi_js = spearman_corr(aphi, js);
    end
    if ~isempty(aphi) && ~isempty(lkl)
        a.spearman_phi_log1pkl = spearman_corr(aphi, lkl);
    end
end

function col = get_col(T, candidates)
    col = [];
    for i = 1:numel(candidates)
        if ismember(candidates{i}, T.Properties.VariableNames)
            col = real(T.(candidates{i}));
            return
        end
    end
end

function m = nanmedian_local(x)
    x = x(isfinite(x));
    if isempty(x), m = nan; else, m = median(x); end
end

function rho = spearman_corr(x, y)
    % Spearman rank correlation without the Statistics Toolbox.
    x = x(:); y = y(:);
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    if numel(x) < 3, rho = nan; return; end
    rx = rank_with_ties(x);
    ry = rank_with_ties(y);
    rx = rx - mean(rx);
    ry = ry - mean(ry);
    denom = sqrt(sum(rx.^2) * sum(ry.^2));
    if denom == 0, rho = nan; else, rho = sum(rx .* ry) / denom; end
end

function r = rank_with_ties(x)
    [~, order] = sort(x);
    n = numel(x);
    r = zeros(n, 1);
    r(order) = 1:n;
    % Average ranks for tied values.
    [xs, ~] = sort(x);
    i = 1;
    while i <= n
        j = i;
        while j < n && xs(j+1) == xs(i)
            j = j + 1;
        end
        if j > i
            avg = mean(i:j);
            r(order(i:j)) = avg;
        end
        i = j + 1;
    end
end

function note_lines = notes_to_strings(notes)
    if iscell(notes)
        note_lines = string(notes(:));
    elseif isstring(notes)
        note_lines = reshape(notes, [], 1);
    elseif ischar(notes)
        note_lines = string(cellstr(notes));
        note_lines = note_lines(:);
    else
        note_lines = reshape(string(notes), [], 1);
    end
end
