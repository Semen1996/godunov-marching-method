function q = flow_speed(u, v)
  % Speed magnitude
  %
  % Input:
  %   u — x velocity component (scalar or array)
  %   v — y velocity component (scalar or array of the same size as u)
  %
  % Output:
  %   q — magnitude of the velocity vector (u, v), q >= 0; array of the same
  %       size as u

  q = sqrt(u.^2 + v.^2);
end
