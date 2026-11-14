function H = total_enthalpy(p, rho, u, v, gamma)
  % Total enthalpy per unit mass of a perfect gas
  %
  % Input:
  %   p     — pressure (scalar or array)
  %   rho   — density (scalar or array of the same size as p)
  %   u     — x velocity component (scalar or array of the same size as p)
  %   v     — y velocity component (scalar or array of the same size as p)
  %   gamma — ratio of specific heats (scalar)
  %
  % Output:
  %   H     — total enthalpy (array of the same size as p)

  H = gamma/(gamma-1)*p./rho + 0.5*(u.^2 + v.^2);
end
