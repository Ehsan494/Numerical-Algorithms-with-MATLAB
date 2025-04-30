% Implicit Backward-Euler
clear; clc;

% Problem parameters
L     = 1.0;            % rod length
alpha = 1e-4;           % thermal diffusivity
Tend  = 500;            % final time

N  = 40;                % spatial sub-intervals
dx = L/N;
dt = 1.0;               % deliberately large to show unconditional stability
Nt = ceil(Tend/dt);     % number of time steps
dt = Tend/Nt;           % adjust to land exactly on Tend
lambda = alpha*dt/dx^2;

% ---- grids & initial/boundary conditions -----------------------------
x  = linspace(0,L,N+1);           % 0 .. L  (N+1 nodes)
%u  = sin(pi*x/L);                 % initial profile  f(x)
u = sin(pi*x/L) + 0.25*sin(6*pi*x/L);
% Dirichlet ends are zero here: u(1)=u(end)=0 for all time

% ---- tridiagonal matrix  A u^{n+1} = u^{n}  --------------------------
main = (1 + 2*lambda) * ones(N-1,1);
off  = -lambda          * ones(N-2,1);
A    = diag(main) + diag(off,1) + diag(off,-1);

% ---- time marching ----------------------------------------------------
for n = 1:Nt
    rhs      = u(2:N)';          % interior old values (column vector)
    u(2:N)   = A \ rhs;          % solve implicit system
    % boundaries remain zero automatically
end

% ---- plot -------------------------------------------------------------
plot(x,u,'LineWidth',2); grid on
xlabel('x'); ylabel('u(x,T_{end})');
title(sprintf('Implicit BTCS (simple):  N=%d,  dt=%.3g,  \\lambda=%.3g',N,dt,lambda));
