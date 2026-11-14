function [WLS, WRS] = riemann_solve(WL, WR, gamma)
  % exact solution of the steady discontinuity breakup problem
  %
  % Lower flow W_L, upper flow W_R. Solution: C- wave in the lower flow, slip
  % line (p*, theta*), C+ wave in the upper flow. The pressure p* is the root of
  %   Theta(p) = f_L(p) + f_R(p) - (theta_L - theta_R) = 0
  % (Newton's method safeguarded by bisection).
  %
  % Input:
  %   WL — [n x 4] lower flow parameters [p, rho, q, theta]: pressure,
  %        density, speed magnitude and its inclination angle to the x axis, rad
  %        (each row is a separate problem, n problems are solved vectorized)
  %   WR — [n x 4] upper flow parameters [p, rho, q, theta]
  %   gamma  — ratio of specific heats (scalar)
  %
  % Output:
  %   WLS — [n x 4] state behind the C- wave (below the slip line),
  %         [p*, rho, q, theta*]
  %   WRS — [n x 4] state behind the C+ wave (above the slip line),
  %         [p*, rho, q, theta*]
  %
  % Errors: shock detachment (required turning exceeds delta_max), vacuum
  % formation, non-positive squared speed behind a wave.


  pL = WL(:,1);
  rhoL = WL(:,2);
  qL = WL(:,3);
  thL = WL(:,4);
  ML = qL./sound_speed(pL, rhoL, gamma);
  HL = gamma/(gamma-1)*pL./rhoL + 0.5*qL.^2;

  pR = WR(:,1);
  rhoR = WR(:,2);
  qR = WR(:,3);
  thR = WR(:,4);
  MR = qR./sound_speed(pR, rhoR, gamma);
  HR = gamma/(gamma-1)*pR./rhoR + 0.5*qR.^2;

  dth = thL - thR;

  n = numel(pL);

  % Upper bound of the search for p*. The pressure behind an attached oblique
  % shock cannot exceed p_k*shock_rmax(M_k): this is the shock with the limiting
  % turning angle delta_max. If p* is larger, the shock detaches. Hence p* is
  % sought on the interval 0 < p* <= hi, where hi is the smaller of the two
  % limiting pressures (section 5.7).
  hi = min(pL.*shock_rmax(ML, gamma), pR.*shock_rmax(MR, gamma));

  if any(wave_f(hi, pL, ML, gamma) + wave_f(hi, pR, MR, gamma) - dth < 0)
    error('riemann_solve: shock detachment (required turning exceeds delta_max)');
  end

  % Vacuum check. wmax = nu(M -> inf) is the limiting Prandtl–Meyer angle
  % through which the flow turns when expanding to p = 0. Flow k can turn
  % a further wmax - nu(M_k), so the flows can diverge by at most
  % 2*wmax - nu(M_L) - nu(M_R). If thR - thL is larger, there is vacuum between them
  % (section 5.7).
  wmax = 0.5*pi*(sqrt((gamma+1)/(gamma-1)) - 1);

  if any(prandtl_meyer_angle(ML, gamma) + prandtl_meyer_angle(MR, gamma) - 2*wmax - dth > 0)
    error('riemann_solve: vacuum formation between the flows');
  end

  % initial guess — linear theory
  lo = zeros(n, 1);
  ZL = sqrt(ML.^2 - 1)./(gamma*pL.*ML.^2);
  ZR = sqrt(MR.^2 - 1)./(gamma*pR.*MR.^2);
  p = (dth + ZL.*pL + ZR.*pR)./(ZL + ZR);

  k = p <= 0;
  p(k) = 0.1*min(pL(k), pR(k));

  k = p >= hi;
  p(k) = 0.5*(lo(k) + hi(k));

  for it = 1:60
    [fL, dfL] = wave_f(p, pL, ML, gamma);
    [fR, dfR] = wave_f(p, pR, MR, gamma);
    T = fL + fR - dth;
    lo(T < 0) = p(T < 0);
    hi(T > 0) = p(T > 0);
    pn = p - T./(dfL + dfR);
    k = ~(pn > lo & pn < hi);
    pn(k) = 0.5*(lo(k) + hi(k));
    dp = max(abs(pn - p)./p);
    p = pn;
    if dp < 1e-12
      break;
    end
  end

  % parameters behind the waves (section 5.5)
  fL = wave_f(p, pL, ML, gamma);
  fR = wave_f(p, pR, MR, gamma);
  thS = 0.5*((thL - fL) + (thR + fR));
  WLS = star_state(p, pL, rhoL, HL, thS, gamma);
  WRS = star_state(p, pR, rhoR, HR, thS, gamma);
end

function W = star_state(pS, pk, rhok, Hk, thS, gamma)
  r = pS./pk;
  rho = rhok.*r.^(1/gamma);                                  % rarefaction
  k = r >= 1;                                                % shock
  rho(k) = rhok(k).*((gamma+1)*r(k) + gamma - 1)./((gamma-1)*r(k) + gamma + 1);
  q2 = 2*Hk - 2*gamma/(gamma-1)*pS./rho;
  if any(q2 <= 0)
    error('riemann_solve: non-positive squared speed behind a wave');
  end
  q = sqrt(q2);
  W = [pS, rho, q, thS];
end
