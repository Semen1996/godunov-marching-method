function c = sound_speed(p, rho, gamma)
  % Speed of sound in a perfect gas
  %
  % Input:
  %   p   — pressure (scalar or array)
  %   rho — density (scalar or array of the same size as p)
  %   gamma   — ratio of specific heats gamma (scalar)
  %
  % Output:
  %   c   — speed of sound, array of the same size as p and rho
  %
  % All quantities in consistent units (dimensionless in main.m).

  c = sqrt(gamma*p./rho);
end
