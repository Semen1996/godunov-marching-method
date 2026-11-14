function M = mach_number(p, rho, u, v, gamma)
  % Mach number
  %
  % Input:
  %   p   — pressure (scalar or array)
  %   rho — density (scalar or array of the same size as p)
  %   u   — x velocity component (scalar or array of the same size as p)
  %   v   — y velocity component (scalar or array of the same size as p)
  %   gamma   — ratio of specific heats (scalar)
  %
  % Output:
  %   M   — Mach number, ratio of the speed magnitude to the speed of sound;
  %         array of the same size as p

  M = flow_speed(u, v)./sound_speed(p, rho, gamma);
end
