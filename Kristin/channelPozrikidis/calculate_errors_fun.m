function [errors] = calculate_errors_fun(u_error,v_error,xx1,xx2,dx,dy,y1,y2,method_name)

scalefactor = 3;
skip = 1;

u_error_max = max(max(abs(u_error)));
v_error_max = max(max(abs(v_error)));
u_error_L1 = dx*dy*sum(sum(abs(u_error)));
v_error_L1 = dx*dy*sum(sum(abs(v_error)));
u_error_L2 = sqrt(dx*dy*sum(sum(u_error.^2)));
v_error_L2 = sqrt(dx*dy*sum(sum(v_error.^2)));
error_mag = sqrt(u_error.^2 + v_error.^2);

ColNames = {'Max','L1','L2'};
Emax = [u_error_max; v_error_max];
EL1 = [u_error_L1; v_error_L1];
EL2 = [u_error_L2; v_error_L2];
errors = table(Emax, EL1, EL2,'VariableNames',ColNames);

figure
subplot(2,2,1);
pcolor(xx1, xx2, u_error)
shading interp
colorbar('southoutside')
title(['Error in u (horizontal) velocity -- ' method_name]);
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
axis equal;
set(gca,'FontSize',12)

subplot(2,2,2);
pcolor(xx1, xx2, v_error)
shading interp
colorbar('southoutside')
title(['Error in v (vertical) velocity -- ' method_name]);
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
axis equal;
set(gca,'FontSize',12)

subplot(2,2,3);
quiver(xx1(1:skip:end,1:skip:end), xx2(1:skip:end,1:skip:end), u_error(1:skip:end,1:skip:end), v_error(1:skip:end,1:skip:end),scalefactor)
hold on
plot(y1,y2,'.k')
axis equal
ylim([-1.1,1.1])
xlim([xx1(1),xx1(end)])
title('Quiver of errors')
set(gca,'FontSize',12)

subplot(2,2,4);
pcolor(xx1, xx2, error_mag)
shading interp
colorbar('southoutside')
title(['Magnitude of Error in Velocity -- ' method_name]);
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
axis equal;
set(gca,'FontSize',12)

end