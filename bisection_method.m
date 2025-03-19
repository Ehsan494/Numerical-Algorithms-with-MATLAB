function [root, error] = bisection_method(f, a, b, tol)
    % Check if the function changes sign
    if f(a) * f(b) >= 0
        error('Invalid interval! The function must change signs.');
    end
    
    error = abs(b - a);
    
    while error >= tol 
        c = (a + b) / 2; % Compute midpoint
        fc = f(c);
        
        % Check if root is found
        if abs(fc) < tol
            root = c;
            return;
        end
        
        % Update the interval
        if f(a) * fc < 0
            b = c;
        else
            a = c;
        end
        
        error = abs(b - a);
        
    end
    
    % Return the final approximate root
    root = (a + b) / 2;
end
