% heat1d_compare_dt.m ---------------------------------------------------
% Goal: examine how the time-step size Δt (dt) simultaneously
%       • affects the stability of the explicit FTCS scheme, and
%       • affects the accuracy of both FTCS (explicit) and BTCS (implicit)
%         for the 1-D heat equation  u_t = α u_xx  with Dirichlet BCs.
% ----------------------------------------------------------------------
clear;  clc;                        % start with a clean workspace & console

%% ---------------------------------------------------------------------
%% 1.  Physical parameters and spatial grid (fixed for all tests)
%% ---------------------------------------------------------------------
L      = 1.0;       % rod length
alpha  = 1e-4;      % thermal diffusivity  (m^2/s)
Tend   = 10.0;      % final simulation time (seconds)

N      = 80;                    % number of spatial sub-intervals
dx     = L/N;                   % spatial step size
x      = linspace(0,L,N+1);     % grid nodes  (row vector)

% Initial temperature profile:
% - Smooth sine wave  +  modest high-frequency bump.
%   The bump seeds a Fourier mode that will blow up quickly if FTCS is unstable.
f = @(x) sin(pi*x/L) + 0.25*sin(6*pi*x/L);

% Analytical solution for the chosen f(x):
% each sine mode decays separately with exp(-α(kπ/L)^2 t).
exact = @(t) ...
         sin(pi*x/L)      .*exp(-alpha*(pi/L   )^2*t) + ...
  0.25 * sin(6*pi*x/L)    .*exp(-alpha*(6*pi/L)^2*t);

dtCrit = 0.5*dx^2/alpha;  % FTCS Courant/Fourier stability limit (λ = 0.5)

%% ---------------------------------------------------------------------
%% 2.  Time-step list to sweep over
%%     (multiples of the critical dt, so λ ranges far below and above 0.5)
%% ---------------------------------------------------------------------
dtList = dtCrit * (0.2 : 0.2 : 50);   % column vector of dt values

%% Arrays that will hold error metrics for each dt
ErrExplicit = NaN(size(dtList));   % max-norm error for FTCS (if stable)
ErrImplicit = NaN(size(dtList));   % max-norm error for BTCS
BlowUp      = false(size(dtList)); % FTCS blew up?  (true/false flag)

%% ---------------------------------------------------------------------
%% 3.  Loop over all candidate dt values
%% ---------------------------------------------------------------------
for k = 1:numel(dtList)
    dt = dtList(k);          % candidate time step chosen from the list
    Nt = ceil(Tend/dt);      % how many whole steps reach or exceed Tend
    % NOTE: we do *not* readjust dt to hit Tend exactly; we want the user-
    %       supplied dt to remain untouched so λ stays as intended.
    
    lambda = alpha*dt/dx^2;  % mesh Fourier number (controls FTCS stability)
    
    % ------------------------------------------------------------------
    % 3a. Explicit FTCS solver for this dt
    % ------------------------------------------------------------------
    u = f(x);            % initialise with the chosen temperature profile
    stable = true;       % assume stable until a NaN/Inf appears
    
    for n = 1:Nt
        unew = u;        % copy old vector (we'll overwrite interior only)
        idx  = 2:N;      % interior node indices (Dirichlet BCs at ends)
        
        % Forward-Time / Centered-Space update
        unew(idx) = u(idx) + lambda * (u(idx+1) - 2*u(idx) + u(idx-1));
        u = unew;        % advance solution
        
        if any(~isfinite(u))      % detect NaN or Inf ⇒ numerical blow-up
            stable = false;
            break                  % exit time loop early
        end
    end
    
    if stable
        ErrExplicit(k) = max(abs(u - exact(Tend)));  % max-norm error
    else
        BlowUp(k) = true;                            % record instability
    end
    
    % ------------------------------------------------------------------
    % 3b. Implicit BTCS solver for the same dt (unconditionally stable)
    % ------------------------------------------------------------------
    u = f(x);                       % reset to initial condition
    
    % Assemble constant tridiagonal matrix A (size [N-1]×[N-1])
    main = (1 + 2*lambda) * ones(N-1,1);
    off  = -lambda         * ones(N-2,1);
    A    = diag(main) + diag(off,1) + diag(off,-1);  % dense for clarity
    
    for n = 1:Nt
        rhs = u(2:N)';           % interior values from previous time level
        u(2:N) = A \ rhs;        % solve A * u^{n+1}_int = u^{n}_int
    end
    
    ErrImplicit(k) = max(abs(u - exact(Tend)));   % max-norm error
end

%% ---------------------------------------------------------------------
%% 4.  Display a table of results
%% ---------------------------------------------------------------------
fprintf('\n%-8s %-10s %-10s %-10s %-8s\n',...
        'dt', 'lambda', 'Err-Expl', 'Err-Impl', 'Blow?');
for k = 1:numel(dtList)
    fprintf('%-8.3g %-10.3f %-10.2e %-10.2e %-8s\n',...
        dtList(k), alpha*dtList(k)/dx^2,...
        ErrExplicit(k), ErrImplicit(k), string(BlowUp(k)));
end

%% ---------------------------------------------------------------------
%% 5.  Log–log plot: error vs dt, highlighting FTCS blow-ups
%% ---------------------------------------------------------------------
figure;
loglog(dtList, ErrImplicit,'-o','LineWidth',1.8); hold on;              % BTCS curve
loglog(dtList(~BlowUp), ErrExplicit(~BlowUp),'-s','LineWidth',1.8);     % FTCS stable curve

% Mark dt values where FTCS blew up with red × on an arbitrary small y-level
if any(BlowUp)
    loglog(dtList(BlowUp), 1e-3*ones(nnz(BlowUp),1),...
           'xr','MarkerSize',10,'LineWidth',1.5);
    legend('Implicit BTCS','Explicit FTCS (stable)','FTCS blew-up');
else
    legend('Implicit BTCS','Explicit FTCS');
end

xlabel('\Delta t');  ylabel('max |u_{num} - u_{exact}|');
title(sprintf('Error vs \\Deltat  (dx = %.3g,  N = %d)', dx, N));
grid on;   set(gca,'XDir','normal');   % ensure left-to-right dt axis
