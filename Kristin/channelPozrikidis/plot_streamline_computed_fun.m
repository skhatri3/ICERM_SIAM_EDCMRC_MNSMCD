function plot_streamline_computed_fun(u1m,u2m,x1gg,x2gg,y1,y2,method_name,beta_value)

plot(y1, y2, 'k.')
hold on

x1gg = x1gg(1:2:end,1:2:end);
x2gg = x2gg(1:2:end,1:2:end);
u1m = u1m(1:2:end,1:2:end);
u2m = u2m(1:2:end,1:2:end);

% Computed streamlines
hh1_comp = streamline(x1gg, x2gg, u1m', u2m', x1gg(1:end, 1), x2gg(1:end, 1)); % streamlines starting at left wall
hh2_comp = streamline(x1gg, x2gg, -u1m', -u2m', x1gg(1:end, end), x2gg(1:end, end)); % streamlines starting at right wall

strmLW = 2;
RGB = orderedcolors("gem");
set(hh1_comp, 'Color', RGB(1,:), 'linewidth', strmLW);
set(hh2_comp, 'Color', RGB(1,:), 'linewidth', strmLW);

title(['Computed streamlines - ' method_name, ' $\beta=$', num2str(beta_value)],'interpreter','latex')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
axis equal
axis([x1gg(1),x1gg(end),x2gg(1)-0.5,x2gg(end)+0.5])
set(gca,'FontSize',12)

end