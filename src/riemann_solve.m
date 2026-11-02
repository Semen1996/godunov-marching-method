function [WLS, WRS] = riemann_solve(WL, WR, gamma)
  % точное решение стационарной задачи о распаде разрыва
  %
  % Нижний поток W_L, верхний W_R. Решение: волна C- в нижнем потоке, линия
  % скольжения (p*, theta*), волна C+ в верхнем. Давление p* — корень
  %   Theta(p) = f_L(p) + f_R(p) - (theta_L - theta_R) = 0
  % (Ньютон с защитой делением пополам).
  %
  % Вход:
  %   WL — [n x 4] параметры нижнего потока [p, rho, q, theta]: давление,
  %        плотность, модуль скорости и угол её наклона к оси x, рад
  %        (каждая строка — отдельная задача, n задач решаются векторно)
  %   WR — [n x 4] параметры верхнего потока [p, rho, q, theta]
  %   gamma  — показатель адиабаты (скаляр)
  %
  % Выход:
  %   WLS — [n x 4] состояние за волной C- (снизу от линии скольжения),
  %         [p*, rho, q, theta*]
  %   WRS — [n x 4] состояние за волной C+ (сверху от линии скольжения),
  %         [p*, rho, q, theta*]
  %
  % Ошибки: отход скачка (требуемый поворот больше delta_max), образование
  % вакуума, неположительный квадрат скорости за волной.


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

  % Верхняя граница поиска p*. Давление за присоединённым косым скачком
  % не может превышать p_k*shock_rmax(M_k): это скачок с предельным углом
  % поворота delta_max. Если p* больше, скачок отходит. Поэтому p* ищем
  % на отрезке 0 < p* <= hi, где hi — меньшее из двух предельных давлений
  % (раздел 5.7).
  hi = min(pL.*shock_rmax(ML, gamma), pR.*shock_rmax(MR, gamma));

  if any(wave_f(hi, pL, ML, gamma) + wave_f(hi, pR, MR, gamma) - dth < 0)
    error('riemann_solve: отход скачка (требуемый поворот больше delta_max)');
  end

  % Проверка на вакуум. wmax = nu(M -> inf) — предельный угол Прандтля–Майера,
  % на который поток поворачивает при расширении до p = 0. Поток k может
  % повернуть ещё на wmax - nu(M_k), значит потоки могут разойтись не более
  % чем на 2*wmax - nu(M_L) - nu(M_R). Если thR - thL больше, между ними вакуум.
  % (раздел 5.7).
  wmax = 0.5*pi*(sqrt((gamma+1)/(gamma-1)) - 1);

  if any(prandtl_meyer_angle(ML, gamma) + prandtl_meyer_angle(MR, gamma) - 2*wmax - dth > 0)
    error('riemann_solve: образование вакуума между потоками');
  end

  % начальное приближение — линейная теория
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

  % параметры за волнами (раздел 5.5)
  fL = wave_f(p, pL, ML, gamma);
  fR = wave_f(p, pR, MR, gamma);
  thS = 0.5*((thL - fL) + (thR + fR));
  WLS = star_state(p, pL, rhoL, HL, thS, gamma);
  WRS = star_state(p, pR, rhoR, HR, thS, gamma);
end

function W = star_state(pS, pk, rhok, Hk, thS, gamma)
  r = pS./pk;
  rho = rhok.*r.^(1/gamma);                                  % разрежение
  k = r >= 1;                                                % скачок
  rho(k) = rhok(k).*((gamma+1)*r(k) + gamma - 1)./((gamma-1)*r(k) + gamma + 1);
  q2 = 2*Hk - 2*gamma/(gamma-1)*pS./rho;
  if any(q2 <= 0)
    error('riemann_solve: неположительный квадрат скорости за волной');
  end
  q = sqrt(q2);
  W = [pS, rho, q, thS];
end
