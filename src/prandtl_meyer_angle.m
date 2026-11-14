function w = prandtl_meyer_angle(M, gamma)
  % Prandtl–Meyer function omega(M) (section 5.3)
  %
  % Input:
  %   M — Mach number (scalar or array, computed element-wise);
  %       for M <= 1 the function returns 0
  %   gamma — ratio of specific heats
  %
  % Output:
  %   w — Prandtl–Meyer angle omega(M) in radians,
  %       array of the same size as M
  k = sqrt((gamma+1)/(gamma-1));
  t = sqrt(max(M.^2 - 1, 0));
  w = k*atan(t/k) - atan(t);
end
