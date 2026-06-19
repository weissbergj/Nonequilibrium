function params = sample_parameter_sets(N, mode)
%SAMPLE_PARAMETER_SETS Log-uniform kinetic rate sampling for a given mode.
%
%   params = sample_parameter_sets(N)          % legacy ranges (default)
%   params = sample_parameter_sets(N, mode)
%
% Columns are [k1, k_1, k2, k_2, k3, k_3]. Ranges come from
% model_constants(mode).rate_log10_bounds:
%   legacy  -> asymmetric per-rate ranges (preserves original behavior)
%   paper_inspired_dimensionless -> symmetric 10^-2..10^2 for all rates
%
% The column-by-column sampling order is preserved so that legacy runs with a
% fixed random seed reproduce the original Overallfinding3 parameter stream.

if nargin < 2 || isempty(mode)
    mode = 'legacy';
end

c = model_constants(mode);
bounds = c.rate_log10_bounds;
if size(bounds, 1) == 1
    bounds = repmat(bounds, 6, 1);
end

params = zeros(N, 6);
for j = 1:6
    lo = bounds(j, 1);
    span = bounds(j, 2) - bounds(j, 1);
    params(:, j) = 10.^(lo + span*rand(N, 1));
end

end
