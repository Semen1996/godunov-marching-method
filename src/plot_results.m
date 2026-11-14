function plot_results(x, ye, yc, Mc, npr, Mj)
  % Mach number field and jet boundary
  %   x  — layer coordinates (n×1)
  %   ye — upper jet boundary on each layer (n×1)
  %   yc — cell centers (n×N)
  %   Mc — Mach number in cells (n×N)
  N = size(Mc, 2);
  XX = repmat(x, 1, N);

  % the jet is long and narrow: wide window, y scale differs from x scale
  figure(1, 'Position', [50, 100, 1400, 550]);  clf;
  contourf(XX, yc, Mc, 30, 'LineColor', 'none');  hold on;
  plot(x, ye, 'k-', 'LineWidth', 1.5);
  hold off;
  colorbar('southoutside');
  axis equal;
  %axis([x(1), x(end), 0, 1.05*max(ye)]);
  set(gca, 'FontSize', 12);
  xlabel('x / h_{noz}');  ylabel('y / h_{noz}');
  title(sprintf('Mach number, M_j = %g, p_j/p_a = %g', Mj, npr));
end
