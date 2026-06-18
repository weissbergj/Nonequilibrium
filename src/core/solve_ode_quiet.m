function [sol, status] = solve_ode_quiet(y0, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3, t_end)
%SOLVE_ODE_QUIET ode15s with per-solve tolerance warning capture.

status = struct('ode_tolerance_warning', false, 'ode_failed', false);
sol = struct('y', []);
lastwarn('');

try
    warn_state = warning('off', 'MATLAB:ode15s:IntegrationTolNotMet');
    c = onCleanup(@() warning(warn_state));
    sol = ode15s(@(t,y) system_ode(t,y,alpha,beta,ATP,k1,k_1,k2,k_2,k3,k_3), ...
        [0 t_end], y0);
catch
    status.ode_failed = true;
    return
end

[msg, warn_id] = lastwarn;
if contains(warn_id, 'IntegrationTolNotMet') || contains(msg, 'Unable to meet integration tolerances')
    status.ode_tolerance_warning = true;
end

if isempty(sol.y) || any(~isfinite(sol.y(:, end)))
    status.ode_failed = true;
end

end
