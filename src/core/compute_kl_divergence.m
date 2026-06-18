function result = compute_kl_divergence(dist_E0, dist_E1, nBins, epsilon, N_stoch)
%COMPUTE_KL_DIVERGENCE Histogram KL between E0 and E1 SSA output samples.
%
% Preserves legacy binning: linspace edges, PDF histcounts, epsilon smoothing.
% Also returns symmetric KL and Jensen-Shannon divergence (supplementary).

if nargin < 5
    N_stoch = numel(dist_E0);
end

result = struct();
result.KL_E0_E1 = nan;
result.KL_E1_E0 = nan;
result.symmetric_KL = nan;
result.JS_divergence = nan;
result.KL_degenerate_histogram = true;
result.KL_sparse_bins = true;
result.insufficient_samples_for_bins = N_stoch < 5*nBins;
result.warning_flags = strings(0,1);

all_data = [dist_E0(:); dist_E1(:)];
all_data = all_data(isfinite(all_data));
if isempty(all_data)
    result.warning_flags(end+1) = "no_finite_distribution_values";
    return
end

degenerate = min(all_data) == max(all_data);
if degenerate
    center = min(all_data);
    half_width = max(epsilon, abs(center)*sqrt(eps));
    edges = linspace(center - half_width, center + half_width, nBins);
else
    edges = linspace(min(all_data), max(all_data), nBins);
end

[p_hist, edges] = histcounts(dist_E0, edges, 'Normalization', 'pdf');
[q_hist, ~] = histcounts(dist_E1, edges, 'Normalization', 'pdf');
bw = edges(2) - edges(1);
p = p_hist*bw + epsilon;
q = q_hist*bw + epsilon;
p = p/sum(p);
q = q/sum(q);

KL_pq = sum(p.*log(p./q));
KL_qp = sum(q.*log(q./p));
m = 0.5*(p + q);
JS = 0.5*sum(p.*log(p./m)) + 0.5*sum(q.*log(q./m));

occupied = sum((p_hist > 0) | (q_hist > 0));
sparse_bins = occupied < max(3, ceil(0.1*(nBins - 1)));

result.KL_E0_E1 = KL_pq;
result.KL_E1_E0 = KL_qp;
result.symmetric_KL = 0.5*(KL_pq + KL_qp);
result.JS_divergence = JS;
result.KL_degenerate_histogram = degenerate;
result.KL_sparse_bins = sparse_bins;

if degenerate
    result.warning_flags(end+1) = "degenerate_histogram";
end
if sparse_bins
    result.warning_flags(end+1) = "sparse_bins";
end
if result.insufficient_samples_for_bins
    result.warning_flags(end+1) = "insufficient_samples_for_bins";
end
if any(~isfinite([KL_pq, KL_qp, JS])) || ~isreal([KL_pq, KL_qp, JS])
    result.warning_flags(end+1) = "nan_or_inf_kl";
end

end
