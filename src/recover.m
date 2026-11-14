function [p, rho, u, v, Dmin] = recover(F, gamma)
  % Conversion from the flux vector F to gas parameters (section 4)
  %
  % Inverse of flux_x. F directly gives v = F3/F1, H = F4/F1 and
  % P = F2/F1 = u + p/(rho*u); for u a quadratic equation remains
  %   (gamma+1)*u^2 - 2*gamma*P*u + 2*(gamma-1)*(H - v^2/2) = 0
  % with discriminant (divided by 4) D = gamma^2*P^2 - 2*(gamma^2-1)*(H - v^2/2).
  % The larger root is taken — it corresponds to a velocity supersonic in x
  % (the smaller one is subsonic). Then p = F2 - F1*u, rho = F1/u.
  %
  % Input:
  %   F    — [n x 4] flux vector through the cross-section x = const
  %          [rho*u, p + rho*u^2, rho*u*v, rho*u*H]
  %          (each row is a separate cell), F(:,1) > 0
  %   gamma — ratio of specific heats (scalar)
  %
  % Output:
  %   p    — [n x 1] pressure
  %   rho  — [n x 1] density
  %   u    — [n x 1] x velocity component
  %   v    — [n x 1] y velocity component
  %   Dmin — smallest value of D over all cells (scalar). If Dmin < 0,
  %          at least one cell has no real root (no gas state corresponds
  %          to such F); the function does not raise an error but computes
  %          u with D = 0 — checking Dmin is left to the calling code.

  v = F(:,3)./F(:,1);
  H = F(:,4)./F(:,1);
  P = F(:,2)./F(:,1);
  D = gamma^2*P.^2 - 2*(gamma^2 - 1)*(H - 0.5*v.^2);
  Dmin = min(D);
  u = (gamma*P + sqrt(max(D, 0)))/(gamma + 1);
  p = F(:,1).*(P - u);
  rho = F(:,1)./u;
end
