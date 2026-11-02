function Phi = face_flux(W, s, gamma, H)
  % поток через грань с наклоном s; H можно задать явно (H_k сохраняется
  % в волнах точно, так поток энергии не зависит от ошибок округления)
  p = W(:,1);  rho = W(:,2);  u = W(:,3);  v = W(:,4);
  if nargin < 4
    H = gamma/(gamma-1)*p./rho + 0.5*(u.^2 + v.^2);
  end
  m = rho.*(v - s.*u);
  Phi = [m, m.*u - s.*p, m.*v + p, m.*H];
end
