function M = isentropic_mach(r, Mk, gamma)
  % Mach number after isentropic expansion to p = r*p_k
  %
  % Input:
  %   r     — pressure ratio p/p_k (scalar or array); r < 1 for expansion
  %   Mk    — Mach number before expansion, at pressure p_k (scalar or array
  %           of the same size as r)
  %   gamma — ratio of specific heats (scalar)
  %
  % Output:
  %   M     — Mach number after expansion to pressure p = r*p_k
  %           (same size as r and Mk)
  %
  % The formula follows from conservation of total temperature and the isentrope:
  %   1 + (gamma-1)/2*M^2 = (1 + (gamma-1)/2*Mk^2) * r^(-(gamma-1)/gamma)

  M = sqrt(2/(gamma-1)*((1 + 0.5*(gamma-1)*Mk.^2).*r.^(-(gamma-1)/gamma) - 1));
end
