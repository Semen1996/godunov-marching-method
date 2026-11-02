function plot_results(x, ye, yc, Mc, npr, Mj)
  % поле числа Маха и граница струи
  %   x  — координаты слоёв (n×1)
  %   ye — верхняя граница струи на каждом слое (n×1)
  %   yc — центры ячеек (n×N)
  %   Mc — число Маха в ячейках (n×N)
  N = size(Mc, 2);
  XX = repmat(x, 1, N);

  % струя длинная и узкая: широкое окно, масштаб по y не равен масштабу по x
  figure(1, 'Position', [50, 100, 1400, 550]);  clf;
  contourf(XX, yc, Mc, 30, 'LineColor', 'none');  hold on;
  plot(x, ye, 'k-', 'LineWidth', 1.5);
  hold off;
  colorbar('southoutside');
  axis equal;
  %axis([x(1), x(end), 0, 1.05*max(ye)]);
  set(gca, 'FontSize', 12);
  xlabel('x / h_{noz}');  ylabel('y / h_{noz}');
  title(sprintf('Число Маха, M_j = %g, p_j/p_a = %g', Mj, npr));
end
