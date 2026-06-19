function E_vals = build_enzyme_sweep(E_min, E_max, n_E_vals, spacing)
%BUILD_ENZYME_SWEEP Total-enzyme sweep vector with selectable spacing.
%
%   E = build_enzyme_sweep(E_min, E_max, n)            % linear (default)
%   E = build_enzyme_sweep(E_min, E_max, n, 'linear')
%   E = build_enzyme_sweep(E_min, E_max, n, 'log')     % logspace
%
% Log spacing uses logspace(log10(E_min), log10(E_max), n) and requires
% E_min > 0. Linear spacing reproduces the original linspace behavior exactly.

if nargin < 4 || isempty(spacing)
    spacing = 'linear';
end
spacing = lower(char(spacing));

switch spacing
    case 'linear'
        E_vals = linspace(E_min, E_max, n_E_vals);
    case 'log'
        if E_min <= 0
            error('build_enzyme_sweep:NonPositiveEmin', ...
                'Log spacing requires E_min > 0 (got %g).', E_min);
        end
        E_vals = logspace(log10(E_min), log10(E_max), n_E_vals);
    otherwise
        error('build_enzyme_sweep:UnknownSpacing', ...
            'Unknown enzyme_spacing "%s". Use "linear" or "log".', spacing);
end

end
