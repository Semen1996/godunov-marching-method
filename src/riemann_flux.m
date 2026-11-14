function Phi = riemann_flux(WL, WR, WLS, WRS, s, gamma)
  % flux through a face from the discontinuity breakup solution (section 5.6)
  %
  % The face is the ray xi = s from the breakup point. The solution is
  % self-similar, the state is constant along a ray; the ray angle selects
  % one of the zones
  %   L | C- wave | L* | slip line | R* | C+ wave | R
  % (a wave is a shock or a rarefaction fan; inside a fan the state is found
  % from the invariant theta -/+ omega(M)), and the flux is computed from the
  % state in that zone. Wave positions are stored as angles (not as slopes
  % tan) so that the comparison with the face ray is unambiguous.
  %
  % Input:
  %   WL  — [n x 4] lower flow parameters [p, rho, q, theta]: pressure,
  %         density, speed magnitude and its inclination angle to the x axis,
  %         rad (each row is a separate face, n faces are computed vectorized);
  %         the flow must be supersonic (M > 1)
  %   WR  — [n x 4] upper flow parameters [p, rho, q, theta], M > 1
  %   WLS — [n x 4] state behind the C- wave (below the slip line),
  %         [p*, rho, q, theta*]
  %   WRS — [n x 4] state behind the C+ wave (above the slip line),
  %         [p*, rho, q, theta*];
  %         [WLS, WRS] = riemann_solve(WL, WR, gamma) for the same WL, WR
  %   s   — [n x 1] face slope dy/dx (tangent of the inclination angle to
  %         the x axis); on a face with s = tan(theta*) the state WLS is taken
  %   gamma — ratio of specific heats (scalar)
  %
  % Output:
  %   Phi — [n x 4] flux through the face per unit length in x,
  %         Phi = G - s*F = [m, m*u - s*p, m*v + p, m*H], m = rho*(v - s*u):
  %         fluxes of mass, x and y momentum and energy (see face_flux).
  %         The total enthalpy H is taken from the original flow on the same
  %         side of the slip line (H_L or H_R) — it is conserved across waves.

  phi = atan(s);
  thL = WL(:,4);  thR = WR(:,4);  thS = WLS(:,4);
  ML = state_mach(WL, gamma);    MR = state_mach(WR, gamma);
  MLS = state_mach(WLS, gamma);  MRS = state_mach(WRS, gamma);
  HL = state_enthalpy(WL, gamma);  HR = state_enthalpy(WR, gamma);

  % wave angles: L1 <= L2 <= theta* <= R1 <= R2 (for a shock L1 = L2, R1 = R2)
  aL1 = thL - asin(1./ML);   aL2 = thS - asin(1./MLS);    % fan: head, tail
  aR1 = thS + asin(1./MRS);  aR2 = thR + asin(1./MR);     % fan: tail, head
  pS = WLS(:,1);  pL = WL(:,1);  pR = WR(:,1);
  k = pS > pL;
  aL1(k) = thL(k) - shock_beta(pS(k)./pL(k), ML(k), gamma);  aL2(k) = aL1(k);
  k = pS > pR;
  aR2(k) = thR(k) + shock_beta(pS(k)./pR(k), MR(k), gamma);  aR1(k) = aR2(k);

  W = WRS;  H = HR;
  k = phi <= thS;  W(k,:) = WLS(k,:);  H(k) = HL(k);
  k = phi <= aL1;  W(k,:) = WL(k,:);
  k = phi >= aR2;  W(k,:) = WR(k,:);
  k = find(phi > aL1 & phi < aL2);
  if ~isempty(k)
    W(k,:) = fan_state(WL(k,:), ML(k), MLS(k), thL(k), HL(k), phi(k), -1, gamma);
  end
  k = find(phi > aR1 & phi < aR2);
  if ~isempty(k)
    W(k,:) = fan_state(WR(k,:), MR(k), MRS(k), thR(k), HR(k), phi(k), +1, gamma);
  end
  q = W(:,3);  th = W(:,4);
  Phi = face_flux([W(:,1), W(:,2), q.*cos(th), q.*sin(th)], s, gamma, H);
end

function M = state_mach(W, gamma)
  % Mach number from the parameter vector
  %
  % Input:
  %   W — [n x 4] flow parameters [p, rho, q, theta]: pressure, density,
  %       speed magnitude and its inclination angle to the x axis, rad
  %   gamma — ratio of specific heats (scalar)
  %
  % Output:
  %   M — [n x 1] Mach number, ratio of the speed magnitude q to the speed of sound

  M = W(:,3)./sound_speed(W(:,1), W(:,2), gamma);
end

function H = state_enthalpy(W, gamma)
  % Total enthalpy per unit mass from the parameter vector
  %
  % Input:
  %   W — [n x 4] flow parameters [p, rho, q, theta]: pressure, density,
  %       speed magnitude and its inclination angle to the x axis, rad
  %   gamma — ratio of specific heats (scalar)
  %
  % Output:
  %   H — [n x 1] total enthalpy H = gamma/(gamma-1)*p/rho + q^2/2

  H = gamma/(gamma-1)*W(:,1)./W(:,2) + 0.5*W(:,3).^2;
end

function b = shock_beta(r, M, gamma)
  % Oblique shock angle to the incoming flow from the pressure ratio
  %
  % Input:
  %   r — pressure ratio across the shock p2/p1, r >= 1 (array)
  %   M — Mach number of the incoming flow, M > 1 (array of the same size as r)
  %   gamma — ratio of specific heats (scalar)
  %
  % Output:
  %   b — angle beta between the shock front and the incoming velocity
  %       vector in radians, from the Mach angle asin(1/M) at r = 1 to pi/2
  %       (normal shock); array of the same size as r:
  %         sin(beta)^2 = (1 + (gamma+1)/(2*gamma)*(r - 1))/M^2
  %       The sine is clipped to one to guard against round-off near a
  %       normal shock.

  b = asin(min(sqrt(1 + (gamma+1)/(2*gamma)*(r - 1))./M, 1));
end

function W = fan_state(Wk, Mk, Mks, thk, Hk, phi, sg, gamma)
  % State inside a rarefaction fan on the ray with angle phi
  %
  % A ray in the fan is a characteristic: phi = theta + sg*mu(M), mu = asin(1/M).
  % Together with the fan invariant theta + sg*omega(M) = theta_k + sg*omega(M_k)
  % this gives an equation for the Mach number on the ray
  %   omega(M) - mu(M) = sg*(theta_k - phi) + omega(M_k),
  % solved by Newton's method safeguarded by bisection on the interval
  % between M_k and M_k*. Then theta = phi - sg*mu(M), pressure and density
  % follow the isentrope from state k, speed — from total enthalpy conservation.
  %
  % Input:
  %   Wk  — [n x 4] flow parameters ahead of the fan [p, rho, q, theta]:
  %         pressure, density, speed magnitude and its inclination angle to
  %         the x axis, rad (each row is a separate ray)
  %   Mk  — [n x 1] Mach number ahead of the fan (at the fan head)
  %   Mks — [n x 1] Mach number behind the fan (at the fan tail); together
  %         with Mk defines the root search interval
  %   thk — [n x 1] velocity inclination angle ahead of the fan, rad (theta from Wk)
  %   Hk  — [n x 1] total enthalpy of the flow ahead of the fan
  %   phi — [n x 1] ray inclination angle to the x axis, rad; must lie between
  %         the fan head and tail
  %   sg  — wave family (scalar): -1 — C- fan (in the lower flow),
  %         +1 — C+ fan (in the upper flow)
  %   gamma — ratio of specific heats (scalar)
  %
  % Output:
  %   W   — [n x 4] state on the ray [p, rho, q, theta]

  rhs = sg*(thk - phi) + prandtl_meyer_angle(Mk, gamma);
  lo = min(Mk, Mks);  hi = max(Mk, Mks);
  M = 0.5*(lo + hi);
  for it = 1:60
    G = prandtl_meyer_angle(M, gamma) - asin(1./M) - rhs;
    lo(G < 0) = M(G < 0);
    hi(G > 0) = M(G > 0);
    t = sqrt(M.^2 - 1);
    dG = t./(M.*(1 + 0.5*(gamma-1)*M.^2)) + 1./(M.*t);
    Mn = M - G./dG;
    k = ~(Mn > lo & Mn < hi);
    Mn(k) = 0.5*(lo(k) + hi(k));
    dM = max(abs(Mn - M)./M);
    M = Mn;
    if dM < 1e-12
      break;
    end
  end
  th = phi - sg*asin(1./M);
  pk = Wk(:,1);  rhok = Wk(:,2);
  p = pk.*((1 + 0.5*(gamma-1)*Mk.^2)./(1 + 0.5*(gamma-1)*M.^2)).^(gamma/(gamma-1));
  rho = rhok.*(p./pk).^(1/gamma);
  q = sqrt(2*Hk - 2*gamma/(gamma-1)*p./rho);
  W = [p, rho, q, th];
end
