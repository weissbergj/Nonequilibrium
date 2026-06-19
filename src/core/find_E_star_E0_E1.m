function out = find_E_star_E0_E1(E_vals, ratio, E0_frac, E1_frac)
%FIND_E_STAR_E0_E1 Midpoint enzyme definition.
%
%   out = find_E_star_E0_E1(E_vals, ratio)                 % legacy 0.9 / 1.1
%   out = find_E_star_E0_E1(E_vals, ratio, E0_frac, E1_frac)
%
% E_star is the swept enzyme value whose output ratio is closest to the
% midpoint of the ratio range. E0 = E0_frac*E_star, E1 = E1_frac*E_star,
% both clipped to the sweep bounds. Defaults preserve legacy behavior
% (E0_frac=0.9, E1_frac=1.1); the dimensionless mode uses 0.95 / 1.05.

if nargin < 3 || isempty(E0_frac)
    E0_frac = 0.9;
end
if nargin < 4 || isempty(E1_frac)
    E1_frac = 1.1;
end

out = struct();
finite_ratio = isfinite(ratio);
if ~any(finite_ratio)
    out.failed = true;
    out.E_star = nan;
    out.E0 = nan;
    out.E1 = nan;
    out.idx_star = nan;
    out.E_star_at_sweep_boundary = false;
    out.small_ratio_dynamic_range = false;
    return
end

y_mid = 0.5*(min(ratio(finite_ratio)) + max(ratio(finite_ratio)));
[~, idx_star] = min(abs(ratio - y_mid));

out.E_star = E_vals(idx_star);
out.E0 = max(E_vals(1), min(E_vals(end), E0_frac * out.E_star));
out.E1 = max(E_vals(1), min(E_vals(end), E1_frac * out.E_star));
out.idx_star = idx_star;
out.failed = false;
out.E_star_at_sweep_boundary = (idx_star == 1) || (idx_star == numel(E_vals));
out.small_ratio_dynamic_range = (max(ratio(finite_ratio)) - min(ratio(finite_ratio))) < 1e-12;

end
