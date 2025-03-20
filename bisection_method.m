function [root, error] = bisection_method(f, a, b, tol, N)

    % Check if the function changes sign
    if f(a) * f(b) >= 0
        error('Invalid interval! The function must change signs.');
    end
    
    error=tol; N=0;
    
    while error >= tol 
        c = (a + b) / 2; % Compute midpoint
        sfc = sign(f(c)); sfa=sign(f(a));
        
        % Update the interval
        if sfa*sfc < 0
            b = c; error = abs(b - a)/2;
        else if sfa*sfc > 0
            a = c; error = abs(b - a)/2;
        else
           root=c; error=0;
        end
        
       N=N+1; 
    end
    
end
