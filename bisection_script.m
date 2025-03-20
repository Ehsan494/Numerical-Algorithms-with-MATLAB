
% A script that finds the approximate root of f(x)=0
% using the bisection method with error less than tol.
%
% Inputs:
%   a, b   - Initial interval [a, b], where f(a) * f(b) < 0
%   tol    - Desired accuracy (stopping criterion)
%   N_max  - Maximum number of iterations to prevent infinite loops
%
% Outputs:
%   root   - Approximated root of f(x) in [a, b]
%   error  - Final error estimate (half the last interval size)
%   iter   - Number of iterations taken 

tol=1e-4; a=-1; b=0; itr=0;
error=(b-a)/2; 

if sign(f(a))* sign(f(b)) >= 0
    error('Invalid interval! The function must change signs over [a, b].');
end 

while error >= tol 
    c=(a+b)/2;
    s=[sign(f(a)) sign(f(c))];

    if prod(s) < 0
        b=c; % root is in [a,c]
        error=(b-a)/2; 
    elseif prod(s) > 0
            a=c; % root is in [c,b]
            error=(b-a)/2; 
    else 
        root=c; error=0; %exact root found
    end
   
  root=c;
  itr=itr+1;
end

fprintf('The results are: root = %.6f, error = %.6f, itrations used = %d\n', root, error, itr);


function [y] = f(x)
y=sin(x)-6*x-5;
end