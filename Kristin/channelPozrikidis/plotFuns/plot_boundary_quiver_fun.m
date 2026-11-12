function plot_boundary_quiver_fun(y1, y2, x1m, x2m, u1_exact, u2_exact, u1m, u2m, skip,scalefactor)

plot(y1,y2,'k.')
hold on
quiver(y1,y2,u1_exact*scalefactor,u2_exact*scalefactor,'r','AutoScale','off')
hold on
quiver(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),u1m(1:skip:end,1:skip:end)*scalefactor,u2m(1:skip:end,1:skip:end)*scalefactor,'k','AutoScale','off')

title('Computed velocities')
%legend('','Exact velocity', 'Computed velocity')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
axis equal
axis([x1m(1),x1m(end)+1,x2m(1)-0.5,x2m(end)+0.5])
set(gca,'FontSize',12)
end