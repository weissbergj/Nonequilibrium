function metrics = compute_steady_state_metrics(Etot, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3, t_end, eps_val, R, T)
%COMPUTE_STEADY_STATE_METRICS Deterministic metrics from one ODE solve.
%
% Formulas match legacy Overallfinding3 local helpers:
% ratio P/(P+S), entropy-like dissipation phi, cycle delta mu, net flux Jnet.

y0 = [0; Etot; 0; 0; 0];
[sol, ode_status] = solve_ode_quiet(y0, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3, t_end);

if ode_status.ode_failed || isempty(sol.y)
    metrics = empty_metrics();
    metrics.ode_tolerance_warning = ode_status.ode_tolerance_warning;
    metrics.ode_failed = ode_status.ode_failed;
    return
end

y = sol.y(:, end);
S = y(1);
E = y(2);
ET = y(3);
ETS = y(4);
P = y(5);

J1f = k1*ATP*E;
J1b = k_1*ET;
J2f = k2*ET*S;
J2b = k_2*ETS;
J3f = k3*ETS;
J3b = k_3*E*P;

flux = [J1f J1b; J2f J2b; J3f J3b];
Phi = 0;
for jj = 1:3
    Jf = flux(jj,1) + eps_val;
    Jb = flux(jj,2) + eps_val;
    Phi = Phi + (Jf - Jb)*R*T*log(Jf/Jb);
end

forward = k1*ATP*k2*S*k3;
backward = k_1*k_2*k_3*P;
dmu = R*T*log((forward + eps_val)/(backward + eps_val));
Jnet = k3*ETS - k_3*E*P;

metrics = struct();
metrics.y = y;
metrics.S = S;
metrics.E = E;
metrics.ET = ET;
metrics.ETS = ETS;
metrics.P = P;
metrics.ratio = P/(P + S + eps_val);
metrics.phi = Phi;
metrics.dmu = dmu;
metrics.Jnet = Jnet;
metrics.flux = flux;
metrics.ode_tolerance_warning = ode_status.ode_tolerance_warning;
metrics.ode_failed = ode_status.ode_failed;

end

function m = empty_metrics()
    m = struct('y', [], 'S', nan, 'E', nan, 'ET', nan, 'ETS', nan, 'P', nan, ...
        'ratio', nan, 'phi', nan, 'dmu', nan, 'Jnet', nan, 'flux', [], ...
        'ode_tolerance_warning', false, 'ode_failed', true);
end
