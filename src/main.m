% MAIN. Planar underexpanded supersonic jet.
%
% A jet with p_j/p_a > 1 exhausts from a planar nozzle of width w_noz into
% a quiescent medium with pressure p_a. The upper half of the jet is computed:
% the symmetry axis below (wall with s_w = 0, section 6.1), the free boundary
% p = p_a above (section 6.3). The Godunov marching method is used with the
% exact solution of the steady discontinuity breakup problem (section 5),
% including on the axis and on the free boundary.

clear;

% --- parameters ---
gamma = 1.4;     % ratio of specific heats
M     = 2.0;     % jet Mach number at the nozzle exit
npr   = 2;       % pressure ratio p_j/p_a
N     = 100;     % number of cells across the half-jet
kappa = 0.8;     % safety factor in the step condition
every = 1;       % print every `every`-th step (0 — do not print)

% geometry (dimensional, m)
w_noz = 0.01;    % half-width of the nozzle exit (from the symmetry plane to the lip)
x_len = 0.2;     % length of the computational domain along the axis

% Dimensionless quantities
y_noz = 1;                 % dimensionless nozzle half-width
x_end = x_len/w_noz;       % dimensionless domain length

pa = 1;                    % dimensionless ambient pressure
pj = npr*pa;
rhoj = 1;
aj = sqrt(gamma*pj/rhoj);  % dimensionless speed of sound
uj = M*aj;

% --- Godunov marching method ---
tic;

p = pj*ones(N, 1);
rho = rhoj*ones(N, 1);
u = uj*ones(N, 1);
v = zeros(N, 1);

eta = linspace(0, 1, N+1)';        % uniform node distribution
x = 0;  ye = y_noz;                % at the exit the upper boundary is the nozzle lip
y = eta*ye;                        % the lower boundary is the axis, y_w = 0


X = [];  YE = [];  MC = [];  YC = [];
n = 0;
while true
  % 1. cell parameters
  c = sound_speed(p, rho, gamma);

  if any(u <= c)
    error('main: u <= c at x = %.4f, marching method not applicable', x);
  end

  % 2. exact breakups on interior faces 3/2 ... N-1/2
  W = [p, rho, flow_speed(u, v), flow_angle(u, v)];   % [p, rho, q, theta]
  WL = W(1:N-1, :);
  WR = W(2:N, :);
  [WLS, WRS] = riemann_solve(WL, WR, gamma);


  % symmetry axis (section 6.1): theta_1 + f_1(p_w) = 0
  pw = axis_boundary(W(1, :), gamma);

  % free boundary (section 6.3): theta_b = theta_N - f_N(p_a)
  thb = free_boundary(W(N, :), pa, gamma);
  sb = tan(thb);

  n = n + 1;
  X(n, 1) = x;  YE(n, 1) = ye;
  MC(n, :) = mach_number(p, rho, u, v, gamma)';  YC(n, :) = 0.5*(y(1:N) + y(2:N+1))';
  if x >= x_end - 1e-12
    break;
  end

  % 3. step (a priori estimate of section 7)
  hc = diff(y);
  q2 = u.^2 + v.^2;
  dx = kappa*min(hc.*(u.^2 - c.^2)./(2*c.*sqrt(q2 - c.^2)));
  dx = min(dx, x_end - x);

  % 4. new layer
  yen = ye + dx*sb;
  yn = eta*yen;
  s = (yn - y)/dx;
  hn = diff(yn);

  % 5–6. fluxes through faces
  Phi = zeros(N+1, 4);
  Phi(1, :) = [0, 0, pw, 0];                   % symmetry axis, s_w = 0
  Phi(2:N, :) = riemann_flux(WL, WR, WLS, WRS, s(2:N), gamma);
  Phi(N+1, :) = [0, -s(N+1)*pa, pa, 0];        % free boundary

  % 7. scheme 3.2
  F = flux_x(p, rho, u, v, gamma);
  Fn = (hc.*F - dx*(Phi(2:N+1, :) - Phi(1:N, :)))./hn;

  % 8. recovery of gas parameters
  [p, rho, u, v, Dmin] = recover(Fn, gamma);
  if Dmin < 0
    error('main: D < 0 at x = %.4f (choking, step too large?)', x + dx);
  end

  x = x + dx;  ye = yen;  y = yn;
  if every > 0 && mod(n, every) == 0
    fprintf('step %6d   x = %9.5f   dx = %.3e   y_e = %.5f\n', n, x, dx, ye);
  end
end

fprintf('Godunov (exact breakup): N = %d, steps %d, time %.1f s\n', N, numel(X) - 1, toc);
plot_results(X, YE, YC, MC, npr, M);
