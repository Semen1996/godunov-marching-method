function th = flow_angle(u, v)
  % Inclination angle of the velocity vector to the x axis
  %
  % Input:
  %   u  — x velocity component (scalar or array)
  %   v  — y velocity component (scalar or array of the same size as u)
  %
  % Output:
  %   th — angle between the vector (u, v) and the x axis in radians, in the
  %        range [-pi, pi]; array of the same size as u and v.
  %        th > 0 for v > 0 (flow deflected upward from the x axis).

  th = atan2(v, u);
end
