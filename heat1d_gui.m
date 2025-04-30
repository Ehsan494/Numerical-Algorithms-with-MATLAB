function heat1d_gui
% heat1d_gui  -----------------------------------------------------------
% Simple dialog-based interface for a 1-D heat-equation demo
% (implicit BTCS, zero-Dirichlet ends)
% ----------------------------------------------------------------------

% ---------- 1. Ask the user for parameters ----------------------------
dlgTitle   = 'Heat-equation parameters';
prompt     = { ...
    'Rod length  L  (m):', ...
    'Thermal diffusivity  \alpha  (m^2/s):', ...
    'Final time  T_end  (s):', ...
    'Spatial cells  N  (integer > 5):', ...
    'Time-step  \Delta t  (s):'};
defVals    = {'1', '1e-4', '5', '40', '1'};   % default values

answer = inputdlg(prompt, dlgTitle, 1, defVals);
if isempty(answer);  return;  end              % user pressed Cancel

% Convert strings to numbers
L      = str2double(answer{1});
alpha  = str2double(answer{2});
Tend   = str2double(answer{3});
N      = round(str2double(answer{4}));
dt     = str2double(answer{5});

% ---------- 2. Sanity checks ------------------------------------------
if any(isnan([L,alpha,Tend,N,dt])) || N<5 || dt<=0 || L<=0 || alpha<=0
    errordlg('Please enter positive numeric values.','Input error');
    return
end

% ---------- 3. Derived grid info --------------------------------------
dx   = L/N;
Nt   = ceil(Tend/dt);          % keep user dt; just pad final step if needed
lambda = alpha*dt/dx^2;

% ---------- 4. Initial & boundary data --------------------------------
x = linspace(0,L,N+1);
u = sin(pi*x/L) + 0.25*sin(6*pi*x/L);  % mixed sine profile

% constant tridiagonal matrix for BTCS
main = (1+2*lambda)*ones(N-1,1);
off  = -lambda      *ones(N-2,1);
A    = diag(main) + diag(off,1) + diag(off,-1);

% ---------- 5. Time stepping (implicit) -------------------------------
for n = 1:Nt
    rhs      = u(2:N)';        % interior old values
    u(2:N)   = A \ rhs;        % solve
end

% ---------- 6. Plot result --------------------------------------------
figure('Name','1-D Heat Equation');
plot(x,u,'LineWidth',2); grid on
xlabel('x (m)'); ylabel('u(x,T_{end})');
title( sprintf(['BTCS result  (N = %d,  dt = %.3g s,  ', ...
                '\\lambda = %.3g)'], N, dt, lambda) );
end
