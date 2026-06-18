function c = model_constants()
%MODEL_CONSTANTS Fixed model parameters (SI units), matching legacy Overallfinding3.

c.alpha = 5.0e-6;
c.beta  = 0.002;
c.ATP   = 1.0e-3;
c.t_end = 10000.0;
c.t_end_stoch_default = 25.0;  % practical default; legacy brute-force used 1500
c.max_steps = 1e6;
c.R = 8.314;
c.T = 300.0;
c.epsilon = 1e-15;
c.E_min = 10e-9;
c.E_max = 5.0e-6;
c.Vol_factor = 1e9;
c.input_signal_definition = 'total_enzyme_concentration';

end
