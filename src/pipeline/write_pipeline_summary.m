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

if exist(screenFile, 'file')
    D = readtable(screenFile);
    summary.screen_total = height(D);
    if ismember('valid_for_screen', D.Properties.VariableNames)
        summary.screen_valid = sum(D.valid_for_screen);
    else
        summary.screen_valid = sum(strcmp(string(D.status), "ok"));
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
end

lines = strings(0, 1);
lines = append_line(lines, "# Pipeline summary");
lines = append_line(lines, "");
lines = append_line(lines, sprintf("- Deterministic rows: **%d**", summary.screen_total));
lines = append_line(lines, sprintf("- Valid for screen: **%d**", summary.screen_valid));
lines = append_line(lines, sprintf("- Selected for SSA: **%d**", summary.selected_total));
lines = append_line(lines, sprintf("- SSA rows: **%d**", summary.ssa_total));
lines = append_line(lines, sprintf("- valid_for_kl: **%d**", summary.ssa_valid_kl));
lines = append_line(lines, sprintf("- max_steps_hit rows: **%d**", summary.ssa_max_steps));
lines = append_line(lines, sprintf("- sparse_bins rows: **%d**", summary.ssa_sparse_bins));
lines = append_line(lines, sprintf("- degenerate_histogram rows: **%d**", summary.ssa_degenerate));
lines = append_line(lines, sprintf("- insufficient_samples_for_bins rows: **%d**", summary.ssa_insufficient_bins));
lines = append_line(lines, sprintf("- failed_trajectories rows: **%d**", summary.ssa_failed_trajectories));
lines = append_line(lines, "");

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
