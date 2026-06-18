function Results = check_metric_refactor_equivalence(N_cases)
%CHECK_METRIC_REFACTOR_EQUIVALENCE Verify shared metrics match legacy formulas.

if nargin < 1 || isempty(N_cases)
    N_cases = 5;
end

mc = model_constants();
Results = table('Size', [N_cases 10], ...
    'VariableTypes', {'double','double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'case_id','phi_abs_diff','dmu_abs_diff','Jnet_abs_diff', ...
    'phi_rel_diff','dmu_rel_diff','Jnet_rel_diff','max_abs_diff','max_rel_diff','passed'});

for case_id = 1:N_cases
    k1  = 10^( 3 + 4*rand);
    k_1 = 10^(-1 + 4*rand);
    k2  = 10^( 3 + 4*rand);
    k_2 = 10^(-1 + 4*rand);
    k3  = 10^(-1 + 4*rand);
    k_3 = 10^( 1 + 4*rand);
    Etot = mc.E_min + (mc.E_max - mc.E_min)*rand;

    legacy_phi = compute_entropy_legacy(Etot, mc, k1, k_1, k2, k_2, k3, k_3);
    legacy_dmu = compute_dmu_legacy(Etot, mc, k1, k_1, k2, k_2, k3, k_3);
    legacy_Jnet = compute_Jnet_legacy(Etot, mc, k1, k_1, k2, k_2, k3, k_3);

    m = compute_steady_state_metrics(Etot, mc.alpha, mc.beta, mc.ATP, ...
        k1, k_1, k2, k_2, k3, k_3, mc.t_end, mc.epsilon, mc.R, mc.T);

    phi_diff = abs(m.phi - legacy_phi);
    dmu_diff = abs(m.dmu - legacy_dmu);
    Jnet_diff = abs(m.Jnet - legacy_Jnet);

    Results.case_id(case_id) = case_id;
    Results.phi_abs_diff(case_id) = phi_diff;
    Results.dmu_abs_diff(case_id) = dmu_diff;
    Results.Jnet_abs_diff(case_id) = Jnet_diff;
    Results.phi_rel_diff(case_id) = rel_diff(m.phi, legacy_phi);
    Results.dmu_rel_diff(case_id) = rel_diff(m.dmu, legacy_dmu);
    Results.Jnet_rel_diff(case_id) = rel_diff(m.Jnet, legacy_Jnet);
    Results.max_abs_diff(case_id) = max([phi_diff, dmu_diff, Jnet_diff]);
    Results.max_rel_diff(case_id) = max([Results.phi_rel_diff(case_id), Results.dmu_rel_diff(case_id), Results.Jnet_rel_diff(case_id)]);
    Results.passed(case_id) = Results.max_abs_diff(case_id) < 1e-10;
end

fprintf('Metric equivalence: %d/%d cases passed (max abs diff < 1e-10).\n', sum(Results.passed), N_cases);
end

function Phi = compute_entropy_legacy(Etot, mc, k1, k_1, k2, k_2, k3, k_3)
    y0 = [0; Etot; 0; 0; 0];
    sol = ode15s(@(t,y) system_ode(t,y,mc.alpha,mc.beta,mc.ATP,k1,k_1,k2,k_2,k3,k_3), [0 mc.t_end], y0);
    y = sol.y(:, end);
    J1f = k1*mc.ATP*y(2); J1b = k_1*y(3);
    J2f = k2*y(3)*y(1); J2b = k_2*y(4);
    J3f = k3*y(4); J3b = k_3*y(2)*y(5);
    Phi = 0;
    for jj = 1:3
        Jf = [J1f J2f J3f]; Jb = [J1b J2b J3b];
        Jf = Jf(jj) + mc.epsilon; Jb = Jb(jj) + mc.epsilon;
        Phi = Phi + (Jf - Jb)*mc.R*mc.T*log(Jf/Jb);
    end
end

function dmu = compute_dmu_legacy(Etot, mc, k1, k_1, k2, k_2, k3, k_3)
    y0 = [0; Etot; 0; 0; 0];
    sol = ode15s(@(t,y) system_ode(t,y,mc.alpha,mc.beta,mc.ATP,k1,k_1,k2,k_2,k3,k_3), [0 mc.t_end], y0);
    y = sol.y(:, end);
    forward = k1*mc.ATP*k2*y(1)*k3;
    backward = k_1*k_2*k_3*y(5);
    dmu = mc.R*mc.T*log((forward + mc.epsilon)/(backward + mc.epsilon));
end

function Jnet = compute_Jnet_legacy(Etot, mc, k1, k_1, k2, k_2, k3, k_3)
    y0 = [0; Etot; 0; 0; 0];
    sol = ode15s(@(t,y) system_ode(t,y,mc.alpha,mc.beta,mc.ATP,k1,k_1,k2,k_2,k3,k_3), [0 mc.t_end], y0);
    y = sol.y(:, end);
    Jnet = k3*y(4) - k_3*y(2)*y(5);
end

function r = rel_diff(a, b)
    denom = max(abs(b), eps);
    r = abs(a - b)/denom;
end
