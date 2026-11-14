function F = flux_x(p, rho, u, v, gamma)
  % Flux vector F through the cross-section x = const (section 1)
  %
  % Input:
  %   p   — [n x 1] pressure (each row is a separate cell)
  %   rho — [n x 1] density
  %   u   — [n x 1] x velocity component
  %   v   — [n x 1] y velocity component
  %   gamma — ratio of specific heats (scalar)
  %
  % Output:
  %   F   — [n x 4] fluxes of mass, x momentum, y momentum and energy through
  %         a unit area of the cross-section x = const:
  %         F = [rho*u, p + rho*u^2, rho*u*v, rho*u*H],
  %         H = gamma/(gamma-1)*p/rho + (u^2 + v^2)/2 — total enthalpy.
  %         The inverse transformation is recover.

  H = gamma/(gamma-1)*p./rho + 0.5*(u.^2 + v.^2);
  m = rho.*u;
  F = [m, p + m.*u, m.*v, m.*H];
end
