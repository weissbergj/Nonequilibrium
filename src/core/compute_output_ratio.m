function [ratio, ode_status] = compute_output_ratio(Etot, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3, t_end, eps_val)
%COMPUTE_OUTPUT_RATIO Final-time P/(P+S) from deterministic ODE.

y0 = [0; Etot; 0; 0; 0];
[sol, ode_status] = solve_ode_quiet(y0, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3, t_end);

if ode_status.ode_failed || isempty(sol.y)
    ratio = nan;
    return
end

y = sol.y(:, end);
ratio = y(5)/(y(5) + y(1) + eps_val);

end
