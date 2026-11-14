function [f, df] = wave_f(p, pk, Mk, gamma)
  % flow turning angle f_k(p) and its derivative (sections 5.2–5.4)
  %
  % Input:
  %   p     — pressure behind the wave (array)
  %   pk    — pressure ahead of the wave, on side k (array of the same size as p)
  %   Mk    — Mach number ahead of the wave, on side k (array of the same size as p)
  %   gamma — ratio of specific heats
  %
  % Output:
  %   f     — flow turning angle f_k(p): oblique shock for p >= pk,
  %           Prandtl–Meyer rarefaction wave for p < pk
  %   df    — derivative df_k/dp

  r = p./pk;

  f = zeros(size(p));
  df = f;

  % oblique shock, p >= p_k
  s = r >= 1;

  M2 = Mk(s).^2;
  rs = r(s);

  a = 2*gamma*M2 - (gamma-1);
  b = gamma + 1;
  c = gamma - 1;

  A  = (rs - 1)./(gamma*M2 - (rs - 1));
  dA = gamma*M2./(gamma*M2 - (rs - 1)).^2;

  B  = max((a - b*rs)./(b*rs + c), 0);
  dB = -b*(a + c)./(b*rs + c).^2;

  sB = sqrt(B);

  t  = A.*sB;
  dt = dA.*sB + A.*dB./(2*max(sB, eps));

  f(s)  = atan(t);
  df(s) = dt./(1 + t.^2)./pk(s);

  % rarefaction wave, p < p_k
  e = ~s;
  M = isentropic_mach(r(e), Mk(e), gamma);
  f(e)  = prandtl_meyer_angle(Mk(e), gamma) - prandtl_meyer_angle(M, gamma);
  df(e) = sqrt(M.^2 - 1)./(gamma*p(e).*M.^2);
end
