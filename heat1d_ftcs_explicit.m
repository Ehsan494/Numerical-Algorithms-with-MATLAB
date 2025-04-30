% heat1d_ftcs_unrestricted.m -------------------------------------------
% Explicit FTCS solver for 1-D heat equation u_t = alpha * u_xx
% with user-selected dt.  Change dt to watch how the solution behaves
% as lambda = alpha*dt/dx^2 crosses the 0.5 stability threshold.
% ----------------------------------------------------------------------
clear; clc;

%% 1. Physical parameters (edit as you like)
L      = 1.0;           % rod length
alpha  = 1e-4;          % thermal diffusivity
Tend   = 500;           % simulation end time (seconds)

%% 2. Grid
N      = 80;            % spatial sub-intervals  (N+1 nodes)
dx     = L / N;

% --- CHOOSE ANY dt YOU WANT HERE --------------------------------------
 dt     = 0.8 * dx^2 / alpha;   % EXAMPLE: lambda = 0.8 (> 0.5) → unstable
%dt   = 0.4 * dx^2 / alpha;   % EXAMPLE: lambda = 0.4 (< 0.5) → stable
 %dt   = 5;                  % another way: huge dt (lambda ≫ 0.5)
% ----------------------------------------------------------------------

Nt     = ceil(Tend / dt);      % number of whole steps
dt     = Tend / Nt;            % tweak dt so Nt*dt = Tend exactly
lambda = alpha * dt / dx^2;
fprintf('\nFTCS run:  N = %d   dt = %.5g   Nt = %d   lambda = %.3f\n',...
        N, dt, Nt, lambda);

%% 3. Initial & boundary conditions
x   = linspace(0, L, N+1);
%u   = sin(pi*x/L);             % initial profile  (smooth; hides instability as dt grows)
% For a worst-case instability test you could try:
%u = (-1).^(0:N);             % checkerboard (+1,-1,+1,-1,...)
% Reasonable initial u
u = sin(pi*x/L) + 0.25*sin(6*pi*x/L);


g0  = @(t) 0;                  % Dirichlet BC at x=0
gL  = @(t) 0;                  % Dirichlet BC at x=L

%% 4. Time stepping loop (explicit FTCS)
for n = 1:Nt
    tnew      = n * dt;

    unew      = u;                       % copy old vector
    unew(1)   = g0(tnew);                % apply BC left
    unew(end) = gL(tnew);                % apply BC right

    idx       = 2:N;                     % interior indices
    unew(idx) = u(idx) + lambda*( u(idx+1) - 2*u(idx) + u(idx-1) );

    u = unew;

    % Optional: stop early if NaNs/Inf appear (blow-up detection)
    if any(~isfinite(u))
        warning('Solution blew up at step %d (t = %.3g s)\n', n, tnew);
        break
    end
end

%% 5. Plot final profile
plot(x, u, 'LineWidth', 2); grid on;
xlabel('x'); ylabel('u(x,T_{end})');
title(sprintf('FTCS   \\lambda = %.3f   (dt = %.3g s)', lambda, dt));
