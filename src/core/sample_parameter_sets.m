function params = sample_parameter_sets(N)
%SAMPLE_PARAMETER_SETS Log-uniform kinetic rate sampling (legacy ranges).

params = zeros(N, 6);
params(:,1) = 10.^( 3 + 4*rand(N,1));
params(:,2) = 10.^(-1 + 4*rand(N,1));
params(:,3) = 10.^( 3 + 4*rand(N,1));
params(:,4) = 10.^(-1 + 4*rand(N,1));
params(:,5) = 10.^(-1 + 4*rand(N,1));
params(:,6) = 10.^( 1 + 4*rand(N,1));

end
