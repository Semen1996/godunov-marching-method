function r = shock_rmax(M, gamma)
  % Limiting pressure ratio on the weak branch of an oblique shock (delta = delta_max)
  %
  % Input:
  %   M     — Mach number of the incoming flow, M > 1 (scalar or array)
  %   gamma — ratio of specific heats (scalar)
  %
  % Output:
  %   r     — pressure ratio p2/p1 across the shock at the maximum flow
  %           turning angle delta_max, array of the same size as M

  M2 = M.^2;
  % sin^2 of the shock angle beta at which the turning angle delta is maximal
  sinB2 = ((gamma+1)*M2 - 4 + sqrt((gamma+1)*((gamma+1)*M2.^2 + 8*(gamma-1)*M2 + 16)))./(4*gamma*M2);
  r = 1 + 2*gamma/(gamma+1)*(M2.*sinB2 - 1);
end
