function c = model_constants(mode)
%MODEL_CONSTANTS Fixed model parameters for a given analysis mode.
%
%   c = model_constants()                          % legacy (default)
%   c = model_constants('legacy')
%   c = model_constants('paper_inspired_dimensionless')
%
% Two modes are supported:
%
%   'legacy'
%       Exactly preserves Ishfaq's original Overallfinding3 values:
%       SI-scale concentrations, ATP=1e-3, R=8.314, T=300, asymmetric
%       log-uniform rate ranges, E0=0.9*E*, E1=1.1*E*.
%
%   'paper_inspired_dimensionless'
%       Same enzyme ODE/SSA model, but dimensionless and paper-inspired:
%       ATP=1, RT=1 (phi and delta_mu reported in kBT-like units),
%       symmetric log-uniform rate ranges 10^-2..10^2, E0=0.95*E*,
%       E1=1.05*E*. This is the mode that avoids the known artifacts.
%
% This does NOT reimplement the PNAS transcriptional Markov model or its
% analytic information-rate derivation. It keeps the enzyme model and only
% changes units, ranges, and signal separation.

if nargin < 1 || isempty(mode)
    mode = 'legacy';
end
mode = normalize_mode(mode);

switch mode
    case 'legacy'
        c.alpha = 5.0e-6;
        c.beta  = 0.002;
        c.ATP   = 1.0e-3;
        c.t_end = 10000.0;
        c.t_end_stoch_default = 25.0;
        c.max_steps = 1e6;
        c.R = 8.314;
        c.T = 300.0;
        c.epsilon = 1e-15;          % ODE flux floor
        c.kl_epsilon = 1e-15;       % histogram KL smoothing
        c.E_min = 10e-9;
        c.E_max = 5.0e-6;
        c.Vol_factor = 1e9;
        c.E0_frac = 0.9;
        c.E1_frac = 1.1;
        % Per-rate log10 bounds [lo hi] for k1,k_1,k2,k_2,k3,k_3.
        c.rate_log10_bounds = [3 7; -1 3; 3 7; -1 3; -1 3; 1 5];
        c.energy_units = 'RT_joule_per_mol';
        c.input_signal_definition = 'total_enzyme_concentration';

    case 'paper_inspired_dimensionless'
        c.alpha = 1.0;
        c.beta  = 1.0;
        c.ATP   = 1.0;
        c.t_end = 10000.0;
        c.t_end_stoch_default = 25.0;
        c.max_steps = 1e6;
        c.R = 1.0;                  % RT = 1 -> phi, delta_mu in kBT-like units
        c.T = 1.0;
        c.epsilon = 1e-15;          % ODE flux floor
        c.kl_epsilon = 1e-15;       % histogram KL smoothing (raw KL kept for continuity)
        c.E_min = 1.0;
        c.E_max = 500.0;
        c.Vol_factor = 1.0;         % E values already molecule-like (1..500)
        c.E0_frac = 0.95;           % delta_c = 0.1*c* (paper convention)
        c.E1_frac = 1.05;
        % Symmetric log-uniform ranges for all six rates: 10^-2 .. 10^2.
        c.rate_log10_bounds = [-2 2];
        c.energy_units = 'kBT';
        c.input_signal_definition = 'total_enzyme_concentration';
        % WARNING: the "input" here is TOTAL ENZYME concentration (E0/E1 around
        % E*). This is NOT the same as the PNAS paper's input concentration c
        % (a signal/ligand level). Varying total enzyme changes the catalyst
        % budget, not an upstream signal, so the discrimination curve is only an
        % ANALOGY to the paper's input-discrimination result.
        %
        % TODO (future, if the researcher confirms the intended signal):
        %   add an input_signal = 'substrate_or_alpha' preset that instead
        %   sweeps substrate influx alpha (or substrate level S) to define
        %   c0/c1, keeping enzyme fixed. That would map more directly onto the
        %   paper's "input concentration" semantics. See docs/MODEL_AUDIT.md.

    otherwise
        error('model_constants:UnknownMode', ...
            'Unknown mode "%s". Use "legacy" or "paper_inspired_dimensionless".', mode);
end

c.mode = mode;

end

function mode = normalize_mode(mode)
    mode = lower(char(mode));
    switch mode
        case {'legacy', 'original'}
            mode = 'legacy';
        case {'paper_inspired_dimensionless', 'dimensionless', 'paper_inspired', 'paper'}
            mode = 'paper_inspired_dimensionless';
    end
end
