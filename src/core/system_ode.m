function dydt = system_ode(~, y, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3)
%SYSTEM_ODE Five-species enzyme cycle ODE (matches legacy Overallfinding3).

S   = y(1);
E   = max(y(2), 0);
ET  = y(3);
ETS = y(4);
P   = y(5);

dS   = alpha - k2*ET*S + k_2*ETS;
dE   = -k1*ATP*E + k_1*ET + k3*ETS - k_3*E*P;
dET  = k1*ATP*E - k_1*ET - k2*ET*S + k_2*ETS;
dETS = k2*ET*S - (k_2 + k3)*ETS + k_3*E*P;
dP   = k3*ETS - k_3*E*P - beta*P;

dydt = [dS; dE; dET; dETS; dP];

end
