f = @(x) sin(x)-6*x-5; % Define function
a = -1; % Lower bound
b = 0; % Upper bound
tol = 1e-4; % Tolerance
N_max = 13; % Maximum iterations

[root, error, iter] = bisection_method(f, a, b, tol, N_max);
fprintf('Root: %.6f, Error: %.6f, Iterations: %d\n', root, error, iter);
