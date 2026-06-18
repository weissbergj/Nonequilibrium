function Results = Overallfinding3(N_runs)
clc; close all;
if nargin==0
    N_runs = 10000;   
end

%% ================= FIXED PARAMETERS (SI UNITS) =================
alpha = 5.0e-6;   
beta  = 0.002;    
ATP   = 1.0e-3;   
t_end = 10000.0;  
t_end_stoch = 1500.0; 
max_steps = 1e6;  

R = 8.314;        
T = 300.0;        
eps_val = 1e-15;  

% Preallocate Results structure
Results = struct('E_vals', cell(1, N_runs), 'ratio', cell(1, N_runs), ...
                 'dist_E0', cell(1, N_runs), 'dist_E1', cell(1, N_runs), ...
                 'k', cell(1, N_runs), 'E_star', cell(1, N_runs), ...
                 'E0', cell(1, N_runs), 'E1', cell(1, N_runs), ...
                 'phi_E0', cell(1, N_runs), 'phi_E1', cell(1, N_runs), ...
                 'dmu_E0', cell(1, N_runs), 'dmu_E1', cell(1, N_runs), ...
                 'KL_pq', cell(1, N_runs), 'KL_qp', cell(1, N_runs), ...
                 'Jnet_E0', cell(1, N_runs), 'Jnet_E1', cell(1, N_runs));

fprintf('Starting 10,000 simulations across 28 parallel CPU workers...\n');

%% =================== MONTE CARLO LOOP (PARFOR) ===================
parfor run = 1:N_runs
    
    % ===== Random rate constants (Biologically scaled log-uniform) =====
    k1  = 10^( 3 + 4*rand );   
    k_1 = 10^( -1 + 4*rand );  
    k2  = 10^( 3 + 4*rand );   
    k_2 = 10^( -1 + 4*rand );  
    k3  = 10^( -1 + 4*rand );   
    k_3 = 10^( 1 + 4*rand );   
    
    %% ================= MIDPOINT DETECT =================
    E_vals = linspace(10e-9, 5.0e-6, 500); 
    ratio  = zeros(size(E_vals));
    
    for i = 1:length(E_vals)
        Etot = E_vals(i);
        y0 = [0; Etot; 0; 0; 0];
        
        sol = ode15s(@(t,y) system_ode(t,y,alpha,beta,ATP,...
            k1,k_1,k2,k_2,k3,k_3), [0 t_end], y0);
        
        % FIX: Safely pull the absolute final state computed by the solver
        y_final = sol.y(:, end); 
        
        S = y_final(1);
        P = y_final(5);
        ratio(i) = P/(P+S+eps_val);
    end
    
    % The ODEs automatically establish the midpoint curve profiles:
    y_mid = 0.5*(min(ratio)+max(ratio));
    [~, idx_star] = min(abs(ratio - y_mid));
    
    E_star = E_vals(idx_star);
    E0 = 0.9 * E_star;
    E1 = 1.1 * E_star;
    
    E0 = max(E_vals(1), min(E_vals(end), E0));
    E1 = max(E_vals(1), min(E_vals(end), E1));
    
    Results(run).E_vals = E_vals;
    Results(run).ratio  = ratio;
    
    %% ================= STOCHASTIC SIMULATIONS =================
    Vol_factor = 1e9; 
    N_stoch = 200; 
    dist_E0_local = zeros(N_stoch,1);
    dist_E1_local = zeros(N_stoch,1);
    
    for i=1:N_stoch
        dist_E0_local(i) = gillespie_local(round(E0 * Vol_factor), Vol_factor, alpha, k1, k2, k_3, k_1, k_2, k3, beta, ATP, t_end_stoch, max_steps, eps_val);
        dist_E1_local(i) = gillespie_local(round(E1 * Vol_factor), Vol_factor, alpha, k1, k2, k_3, k_1, k_2, k3, beta, ATP, t_end_stoch, max_steps, eps_val);
    end
    Results(run).dist_E0 = dist_E0_local;
    Results(run).dist_E1 = dist_E1_local;
 
    %% ================= ENERGY CALCULATIONS =================
    Results(run).phi_E0 = compute_entropy_local(E0, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3, t_end, eps_val, R, T);
    Results(run).phi_E1 = compute_entropy_local(E1, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3, t_end, eps_val, R, T);
    
    %% ================= CHEMICAL POTENTIAL =================
    Results(run).dmu_E0 = compute_dmu_local(E0, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3, t_end, eps_val, R, T);
    Results(run).dmu_E1 = compute_dmu_local(E1, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3, t_end, eps_val, R, T);
    
    %% ================= KL DIVERGENCE =================
    all_data = [dist_E0_local; dist_E1_local];
    bins = linspace(min(all_data), max(all_data), 40);
    [p_hist, edges] = histcounts(dist_E0_local, bins, 'Normalization', 'pdf');
    [q_hist, ~]     = histcounts(dist_E1_local, bins, 'Normalization', 'pdf');
    
    bw = edges(2)-edges(1);
    p = p_hist*bw + eps_val; 
    q = q_hist*bw + eps_val;
    p = p/sum(p); 
    q = q/sum(q);
    
    Results(run).KL_pq = sum(p.*log(p./q));
    Results(run).KL_qp = sum(q.*log(q./p));
    
    %% ================= VERIFY NET FLUXES =================
    Results(run).Jnet_E0 = compute_Jnet_local(E0, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3, t_end);
    Results(run).Jnet_E1 = compute_Jnet_local(E1, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3, t_end);
    
    %% ================= STORE CONSTANTS =================
    Results(run).k = [k1 k_1 k2 k_2 k3 k_3];
    Results(run).E_star = E_star;
    Results(run).E0 = E0;
    Results(run).E1 = E1;
    
    if mod(run, 50) == 0
        fprintf('Worker thread registered milestone at iteration %d\n', run);
    end
end

save('MonteCarloResults.mat', 'Results')
fprintf('\n==== Monte Carlo Complete! Saved all 10000 runs. ====\n')
end

%% ====================================================
%% =========== PARFOR COMPATIBLE HELPERS ==============
%% ====================================================
function Phi = compute_entropy_local(Etot, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3, t_end, eps_val, R, T)
    y0=[0;Etot;0;0;0];
    sol=ode15s(@(t,y) system_ode(t,y,alpha,beta,ATP,...
        k1,k_1,k2,k_2,k3,k_3),[0 t_end],y0);
    y=sol.y(:, end); % FIX Applied
    S=y(1); E=y(2); ET=y(3); ETS=y(4); P=y(5);
    J1f=k1*ATP*E; J1b=k_1*ET;
    J2f=k2*ET*S;  J2b=k_2*ETS;
    J3f=k3*ETS;   J3b=k_3*E*P;
    flux=[J1f J1b; J2f J2b; J3f J3b];
    Phi=0;
    for jj=1:3
        Jf=flux(jj,1)+eps_val;
        Jb=flux(jj,2)+eps_val;
        Phi=Phi+(Jf-Jb)*R*T*log(Jf/Jb);
    end
end

function dmu = compute_dmu_local(Etot, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3, t_end, eps_val, R, T)
    y0=[0;Etot;0;0;0];
    sol=ode15s(@(t,y) system_ode(t,y,alpha,beta,ATP,...
        k1,k_1,k2,k_2,k3,k_3),[0 t_end],y0);
    y=sol.y(:, end); % FIX Applied
    S=y(1); P=y(5);
    forward=k1*ATP*k2*S*k3;
    backward=k_1*k_2*k_3*P;
    dmu=R*T*log((forward+eps_val)/(backward+eps_val));
end

function Jnet = compute_Jnet_local(Etot, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3, t_end)
    y0=[0;Etot;0;0;0];
    sol=ode15s(@(t,y) system_ode(t,y,alpha,beta,ATP,...
        k1,k_1,k2,k_2,k3,k_3),[0 t_end],y0);
    y=sol.y(:, end); % FIX Applied
    E=y(2); ETS=y(4); P=y(5);
    Jnet=k3*ETS-k_3*E*P;
end

function value = gillespie_local(Etot_molecules, Vol, alpha, k1, k2, k_3, k_1, k_2, k3, beta, ATP, t_end_stoch, max_steps, eps_val)
    alpha_mol = alpha * Vol;
    k1_vol = k1 / Vol;
    k2_vol = k2 / Vol;
    k_3_vol = k_3 / Vol;
    ATP_mol = round(ATP * Vol);
    
    S=0; E=Etot_molecules; ET=0; ETS=0; P=0; t=0;
    steps = 0; 
    
    while (t < t_end_stoch) && (steps < max_steps)
        steps = steps + 1;
        
        a=[alpha_mol;
           k1_vol*ATP_mol*E;
           k_1*ET;
           k2_vol*ET*S;
           k_2*ETS;
           k3*ETS;
           k_3_vol*E*P;
           beta*P];
        a0=sum(a);
        if a0<=0, break, end
        
        tau=-log(rand)/a0;
        t=t+tau;
        r=find(cumsum(a)>=rand*a0,1);
        
        switch r
            case 1, S=S+1;
            case 2, if E>0,   E=E-1;   ET=ET+1; end
            case 3, if ET>0,  ET=ET-1;  E=E+1;   end
            case 4, if ET>0 && S>0, ET=ET-1; S=S-1; ETS=ETS+1; end
            case 5, if ETS>0, ETS=ETS-1; ET=ET+1; S=S+1; end
            case 6, if ETS>0, ETS=ETS-1; E=E+1;   P=P+1; end
            case 7, if E>0 && P>0,   E=E-1;   P=P-1;   ETS=ETS+1; end
            case 8, if P>0,   P=P-1;   end
        end
    end
    value=P/(P+S+eps_val);
end

%% ================= EXTERNAL ODE SYSTEM =================
function dydt = system_ode(~, y, alpha, beta, ATP, k1, k_1, k2, k_2, k3, k_3)
S   = y(1);
E   = max(y(2),0);
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
