function thb = free_boundary(WN, pa, gamma)
  % free boundary (section 6.3): C- wave in cell N, p* = p_a
  %
  % Input:
  %   WN — parameter vector in cell N, adjacent to the jet boundary:
  %        [p, rho, q, theta] — pressure, density, speed magnitude and its
  %        inclination angle to the x axis, rad; the flow must be supersonic (M > 1)
  %   pa — ambient pressure (pressure behind the C- wave at the boundary)
  %   gamma — ratio of specific heats (scalar)
  %
  % Output:
  %   thb — inclination angle of the jet boundary to the x axis in radians —
  %         the flow angle behind the C- wave (a shock for pa > p_N, a
  %         rarefaction wave for pa < p_N):
  %         thb = theta_N - f_N(pa)
  %
  % Error if pa exceeds the limiting pressure behind an attached oblique
  % shock (p_N * shock_rmax(M_N, gamma)) — shock detachment.
  thN = WN(4);
  MN = WN(3)/sound_speed(WN(1), WN(2), gamma);
  if pa > WN(1)*shock_rmax(MN, gamma)
    error('free_boundary: shock detachment at the jet boundary');
  end
  thb = thN - wave_f(pa, WN(1), MN, gamma);
end
