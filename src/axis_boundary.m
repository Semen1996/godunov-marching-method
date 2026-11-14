function pw = axis_boundary(W1, gamma)
  % Pressure on the symmetry axis (section 6.1): wall with s_w = 0, C+ wave in cell 1.
  % Behind the wave the flow runs along the axis, theta* = 0:
  %   theta_1 + f_1(p_w) = 0
  % Reduces to a discontinuity breakup with cell 1 mirrored about the axis:
  % W_L = [p1, rho1, q1, -theta1], W_R = W1. Then theta_L = -theta_1,
  % f_L = f_R = f_1, and the riemann_solve equation
  %   f_L + f_R - (theta_L - theta_R) = 2*(f_1 + theta_1) = 0
  % coincides with the axis condition, while theta* = 0 by symmetry.
  %
  % Input:
  %   W1    — parameter vector in cell 1, adjacent to the axis:
  %           [p, rho, q, theta] — pressure, density, speed magnitude and
  %           its inclination angle to the x axis, rad;
  %           the flow must be supersonic (M > 1)
  %   gamma — ratio of specific heats (scalar)
  %
  % Output:
  %   pw    — pressure on the axis behind the C+ wave (shock or rarefaction
  %           wave) at which the flow turns along the axis (theta* = 0)
  %
  % Errors (from riemann_solve):
  %   shock detachment from the axis — required turning exceeds delta_max;
  %   vacuum formation — required turning in the rarefaction wave is too large.

  W0 = [W1(1), W1(2), W1(3), -W1(4)];   % mirror cell below the axis
  [~, WRS] = riemann_solve(W0, W1, gamma);
  pw = WRS(1);
end
