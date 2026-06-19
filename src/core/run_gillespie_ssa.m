function [value, steps, hit_max_steps, failed, P_final, S_final] = run_gillespie_ssa(Etot_molecules, Vol, alpha, k1, k2, k_3, k_1, k_2, k3, beta, ATP, t_end_stoch, max_steps, eps_val)
%RUN_GILLESPIE_SSA Gillespie SSA for enzyme model (legacy Overallfinding3).
%
% Extra outputs P_final, S_final return the raw final product and substrate
% molecule counts (for diagnosing whether the P/(P+S) ratio itself causes
% boundary pile-up). Existing callers that request only the first four outputs
% are unaffected.

value = nan;
steps = 0;
hit_max_steps = false;
failed = false;
P_final = nan;
S_final = nan;

try
    alpha_mol = alpha * Vol;
    k1_vol = k1 / Vol;
    k2_vol = k2 / Vol;
    k_3_vol = k_3 / Vol;
    ATP_mol = round(ATP * Vol);

    S = 0; E = Etot_molecules; ET = 0; ETS = 0; P = 0; t = 0;

    while (t < t_end_stoch) && (steps < max_steps)
        steps = steps + 1;
        a = [alpha_mol;
             k1_vol*ATP_mol*E;
             k_1*ET;
             k2_vol*ET*S;
             k_2*ETS;
             k3*ETS;
             k_3_vol*E*P;
             beta*P];
        a0 = sum(a);
        if a0 <= 0
            break
        end

        tau = -log(rand)/a0;
        t = t + tau;
        r = find(cumsum(a) >= rand*a0, 1);

        switch r
            case 1, S = S + 1;
            case 2, if E > 0, E = E - 1; ET = ET + 1; end
            case 3, if ET > 0, ET = ET - 1; E = E + 1; end
            case 4, if ET > 0 && S > 0, ET = ET - 1; S = S - 1; ETS = ETS + 1; end
            case 5, if ETS > 0, ETS = ETS - 1; ET = ET + 1; S = S + 1; end
            case 6, if ETS > 0, ETS = ETS - 1; E = E + 1; P = P + 1; end
            case 7, if E > 0 && P > 0, E = E - 1; P = P - 1; ETS = ETS + 1; end
            case 8, if P > 0, P = P - 1; end
        end
    end

    hit_max_steps = steps >= max_steps;
    value = P/(P + S + eps_val);
    P_final = P;
    S_final = S;
catch
    failed = true;
end

end
