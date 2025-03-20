function [root, error] = bisection_method(f, a, b, tol, N)

    % Check if the function changes sign
    if f(a) * f(b) >= 0
        error('Invalid interval! The function must change signs.');
    end
    
    error=tol; N=0;
    
    while error >= tol 
        c = (a + b) / 2; % Compute midpoint
        fc = f(c);
        
        % Update the interval
        if f(a) * fc < 0
            b = c;
        else if f(b)*fc < 0
            a = c;
        else
           root=c; error=0;
           break; 
        end
        
        error = abs(b - a)/2;
       N=N+1; 
    end
    
end
