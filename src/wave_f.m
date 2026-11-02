function [f, df] = wave_f(p, pk, Mk, gamma)
  % угол поворота f_k(p) и его производная (разделы 5.2–5.4)
  %
  % Входные аргументы:
  %   p     – давление за волной (массив)
  %   pk    – давление перед волной, со стороны k (массив того же размера, что p)
  %   Mk    – число Маха перед волной, со стороны k (массив того же размера, что p)
  %   gamma – показатель адиабаты
  %
  % Выходные аргументы:
  %   f     – угол поворота потока f_k(p): косой скачок при p >= pk,
  %           волна разрежения Прандтля–Майера при p < pk
  %   df    – производная df_k/dp

  r = p./pk;

  f = zeros(size(p));
  df = f;

  % косой скачок, p >= p_k
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

  % волна разрежения, p < p_k
  e = ~s;
  M = isentropic_mach(r(e), Mk(e), gamma);
  f(e)  = prandtl_meyer_angle(Mk(e), gamma) - prandtl_meyer_angle(M, gamma);
  df(e) = sqrt(M.^2 - 1)./(gamma*p(e).*M.^2);
end
