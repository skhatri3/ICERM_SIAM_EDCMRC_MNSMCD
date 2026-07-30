function plot_velocity_pcolor_fun(y1,y2,x1m,x2m,u1m,u2m,ummag,skip,scalefactor)

plot(y1,y2,'k.')
hold on
pc = pcolor(x1m,x2m,ummag);
set(pc,'facealpha',0.8)
shading interp
colorbar('southoutside')
set(gca,'FontSize',12)

axis equal
quiver(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),u1m(1:skip:end,1:skip:end),u2m(1:skip:end,1:skip:end),scalefactor,'k','AutoScale','off')
colorbar
title('Computed Velocity')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
set(gca,'FontSize',12)
end