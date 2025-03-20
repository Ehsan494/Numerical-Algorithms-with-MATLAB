function [root, error, iter] = bisection_method(f, a, b, tol, N_max)
% A MATLAB function to compute an approximate root of f(x) = 0
% using the bisection method with error less than tol.
%
% Inputs:
%   f      - Function handle for f(x)
%   a, b   - Initial interval [a, b], where f(a) * f(b) < 0
%   tol    - Desired accuracy (stopping criterion)
%   N_max  - Maximum number of iterations to prevent infinite loops
%
% Outputs:
%   root   - Approximated root of f(x) in [a, b]
%   error  - Final error estimate (half the last interval size)
%   iter   - Number of iterations taken

    % Check if the function changes sign in the interval
    if f(a) * f(b) >= 0
        error('Invalid interval! The function must change signs over [a, b].');
    end

    % Initialize iteration counter and error estimate
    iter = 0;
    error = abs(b - a) / 2;

    % Start the bisection loop
    while error >= tol && iter < N_max
        c = (a + b) / 2; % Compute midpoint
        fc = f(c);
        
        s=[sign(f(a)) sign(f(c))];

        % Update the interval based on function signs
        if s(1) * s(2) < 0
            b = c; % Root is in [a, c]
        else if s(1)*s(2) > 0
            a = c; % Root is in [c, b]
        else
            root=c; error=0; return; % stop immediately if root is found    
        end

        % Update error estimate and iteration counter
        error = abs(b - a) / 2;
        iter = iter + 1;
    end

    % Final root approximation
    root = c;

    % If max iterations reached, issue a warning
    if iter >= N_max
        warning('Maximum number of iterations reached. Solution may not be accurate.');
    end
end
