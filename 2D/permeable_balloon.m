%main_example3

%Based on Example 3 Section 4.2 of Cortez, SIAM J. Sci Comput. 2001
%Modified extensively
%We have a time dependent problem where the boundary is initialized as an
%ellipse
%Setting forces on a cylinder of radius 1 and computing the
%velocity at given points

%Developed for ICERM
%July 2024
%addpath('/Users/wanda/Documents/GitHub/Regularized_Stokeslets/2D');
%rmpath('/Users/wanda/Documents/GitHub/Regularized_Stokeslets/2D');

close all
clear all
clc

%% Plotting figures
skip = 4; %for quiver plots
set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',2.0,...
    'defaultlinelinewidth',2.0,'defaultlinemarkersize',10.0)

video_flag = 0; % save the movie
if(video_flag)
    vidObj = VideoWriter('balloon_problem.mp4','MPEG-4');
    open(vidObj);
end


%% Parameters to set

%setting the viscosity
mu = 1;

%setting major and minor axis of the initial ellipse
a = 2;
b = 2;
r = sqrt(a*b);

blob_num = 2;

%domain on which velocity is computed
x1min = -4;
x1max = 4;
%x2min and x2max to compare with plots in paper
%x2min = 3/10;
%x2max = 3/10;
%x2min and x2max to see a 2d domain
x2min = -4;
x2max = 4;

%number of points on boundary where force is applied
N = 100;

%resolution for velocity grid
Nx1 = 80;
%Nx2 = 1; %corresponds to compare with plots in paper above
Nx2 = 80; %for a 2d domain

%discretization of cylinder boundary
ds = 2*pi*r/N;
s = 0:ds:2*pi*r;
s = s(1:end-1)';

slopearray = zeros(6,1);
% Loop over epsilon size
%for j=1:6


    %regularization parameter
    ep = 2*ds;
    alpha = 1e-1; % was 1e-2
    %beta = 0.5*ep*alpha*ones(N,1); % check this
    beta = 4/3*ep*alpha*ones(N,1); % check this

    %initial position cylinder on which forces are applied - cylinder is of radius r
    y1 = a*cos(s/r);
    y2 = b*sin(s/r);

    y10 = a*cos(s/r);
    y20 = b*sin(s/r);

    %parameters for tension T = gamma + kstiff(|x_s|-1)
    gamma = 1;
    kstiff = 1e-4*0;

    %N by N difference matrix to compute d/ds        %%%%force f = d\ds(T*tangent)
    e1 = ones(N,1);
    Dc = spdiags([-e1 e1],[0  1],N,N);
    Dc(N,1) =  1;
    Dc = Dc./ds;

    %reference configuration of cylinder
    %resting radius of resulting circle
    y1o = r*cos(s/r);
    y2o = r*sin(s/r);

    yo = [y1o, y2o];
    dyo = Dc*yo;
    magdyo = sqrt(dyo(:,1).^2 + dyo(:,2).^2);

    % figure(1)
    % plot(y1,y2,'k.-')
    % hold on
    % plot(y1o,y2o,'r.-')
    % axis equal

    %figure(2)
    %plot(magdyo)

    %points on which velocity will be computed
    xx1 = linspace(x1min,x1max,Nx1);
    xx2 = linspace(x2min,x2max,Nx2);
    [x1m,x2m] = ndgrid(xx1,xx2);
    x1 = x1m(:);
    x2 = x2m(:);

    %final time
    tfinal = 5;

    %number of timesteps
    tsteps = 500*0+100;

    %timestep
    dtstep = tfinal/tsteps;
    counter = 1;
    area = zeros(tsteps+1,1);
    area(1) = polyarea([y1;y1(1)],[y2;y2(1)]);
    time = zeros(tsteps+1,1);
    time(1) = 0;

    %plot(y1,y2,'k.-')
    %looping over time
    for t = 0:dtstep:tfinal-dtstep

        y = [y1, y2];

        %compute tension
        dy = Dc*y;
        magdy = sqrt(dy(:,1).^2 + dy(:,2).^2);
        ten = gamma + kstiff*(magdy./magdyo - 1);
        dym = -Dc'*y;
        magdym = sqrt(dym(:,1).^2 + dym(:,2).^2);
        maddyom = [magdyo(end);magdyo(1:end-1)];
        tenm = gamma + kstiff*(magdym./maddyom - 1);
        f = (ten.*dy - tenm.*dym)./ds;
        f = f*ds;
        magf = sqrt(f(:,1).^2 + f(:,2).^2);

        tangent = 0.5*(dym+dy);
        magtangent = sqrt(tangent(:,1).^2 + tangent(:,2).^2);
        tangent = tangent./magtangent;
        normal = [tangent(:,2), -tangent(:,1)];
        % quiver(y1,y2,normal(:,1),normal(:,2))
        % axis equal

        ub = RegStokeslets2D_forcetovelocity(y,f,y,ep,mu);
        magub = sqrt(ub(:,1).^2 + ub(:,2).^2);

        u_beta = RegStokeslets2D_permeable_gtovelocity(y,f,...
            y,ep,mu,blob_num, beta, normal);
        magu_beta = sqrt(u_beta(:,1).^2 + u_beta(:,2).^2);

        y1 = y1 + dtstep*(ub(:,1)-u_beta(:,1));
        y2 = y2 + dtstep*(ub(:,2)-u_beta(:,2));

        %computing velocity
        u = RegStokeslets2D_forcetovelocity(y,f,[x1,x2],ep,mu);
        u1 = u(:,1);
        u2 = u(:,2);
        u1m = reshape(u1,size(xx1,2),size(xx2,2));
        u2m = reshape(u2,size(xx1,2),size(xx2,2));

        % compute the pressure
        % [p] = RegStokeslets2D_forcetopressure(y,f,[x1,x2],ep,blob_num);
        % pmesh = reshape(p,size(xx1,2),size(xx2,2));
        % jumpp = max(p)-min(p)
 
        % figure(1)
        % plot(y1,y2,'k.-')
        % hold on
        % plot(y1o,y2o,'r.-')
        % axis equal
        % quiver(y1,y2,f(:,1),f(:,2))

        figure(1)
        hold off
        plot(y1,y2,'k.-');
        hold on
        plot(y10,y20,'r--');
        axis equal
        quiver(y1,y2,u_beta(:,1),u_beta(:,2))
        quiver(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),u1m(1:skip:end,1:skip:end),u2m(1:skip:end,1:skip:end),'k')
        xlim([x1min,x1max])
        ylim([x2min,x2max])
        time_str = sprintf('Time %.1f\n',time(counter));
        title(time_str);

        % figure(1)
        % pfig = pcolor(x1m,x2m,pmesh);
        % shading interp;
        % set(pfig, 'EdgeColor', 'none');
        % hold on
        % plot(y1,y2,'k.-')
        % colorbar
        % time_str = sprintf('Time %.1f\n',time(counter));
        % title(time_str);
        % axis equal
        % 
        % figure(2)
        % plot(x2m(40,:),pmesh(40,:));
        % grid on

        % if(video_flag)
        %     currFrame = getframe(gcf);
        %     writeVideo(vidObj,currFrame);
        % end

        counter = counter + 1;

        area(counter) = polyarea([y1;y1(1)],[y2;y2(1)]);
        time(counter) = time(counter-1) + dtstep;


    end

    if(video_flag)
        close(vidObj);
    end



    %figure;
    %plot(time,area)

    % figure;
    % plot(time,(area-area(1))/area(1))
    % grid on
    % hold on
    % 
    % ind = 101;
    % delta_area = (area-area(1))/area(1);
    % pfit = polyfit(time(ind:end),delta_area(ind:end),1);
    % py = polyval(pfit,time);
    % plot(time,py,'--');
    % hold off
    % slopearray(j) = pfit(1);
    % legend('Displacement','Linear fit');
    % slope = sprintf('Slope =  %g',p(1));
    % text(1,-1e-3,slope,'FontSize',20);
    % pause(0.1);

% end
% 
% 
% figure;
% plot(slopearray,'ok-','markersize',10,'markersize',10,'linewidth',1,'markerfacecolor','#A7A5A5');
% grid on
% xlabel('\epsilon \times \Delta s');
% ylabel('Leakage rate');
%   set(gcf,'color','w');
% saveas(gcf, 'area_loss_rate_summary_newblob', 'epsc');
% 
% %
% %tangential of cylinder boundary
% yp1 = -sin(t);
% yp2 = cos(t);
%
% %forces on cylinder boundary
% %note that the force density is given in the paper - multiply by radius*dt
% f1 = 2*sin(3*t).*yp1*dt;
% f2 = 2*sin(3*t).*yp2*dt;
%

%
% %computing velocity
% u = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[x1,x2],ep,mu);
% u1 = u(:,1);
% u2 = u(:,2);
% u1m = reshape(u1,size(xx1,2),size(xx2,2));
% u2m = reshape(u2,size(xx1,2),size(xx2,2));
%
% %% Computing error
%
% %exact solution
%
% for i = 1:length(xx1)
%
%     for j = 1:length(xx2)
%
%         r = sqrt(x1m(i,j).^2 + x2m(i,j).^2); %radius
%         s = atan2(x2m(i,j),x1m(i,j)); %angle
%
%         if (r < 1)
%
%             uexact1(i,j) = cos(2*s)*r^2/8 + cos(4*s)*r^4/16 - cos(2*s)*r^4/4;
%             uexact2(i,j) = -sin(2*s)*r^2/8 + sin(4*s)*r^4/16 + sin(2*s)*r^4/4;
%
%         else
%
%            uexact1(i,j) = -cos(2*s)/(r^2)/8 + 5*cos(4*s)/(r^4)/16 - cos(4*s)/(r^2)/4;
%            uexact2(i,j) = sin(2*s)/(r^2)/8 + 5*sin(4*s)/(r^4)/16 - sin(4*s)/(r^2)/4;
%
%         end
%
%     end
%
% end
%
% %computing error
% error1 = abs(u1m-uexact1);
% error2 = abs(u2m-uexact2);
% errormag = sqrt(error1.^2 + error2.^2);
%
% %prints the max error
% fprintf('maximum error in u1: %d \n',max(max(error1)));
% fprintf('maximum error in u2: %d \n',max(max(error2)));
%

%
% %plots as in paper, velocity is computed on a line, x2 = constant
% if (Nx2 == 1)
%
%     figure(1)
%     plot(x1,u1,'k.-')
%     hold on
%     plot(x1,uexact1,'r--')
%     plot(x1,u2,'b.-')
%     plot(x1,uexact2,'g--')
%     title('Numerical and Exact Solutions')
%     legend('u1','uexact1','u2','uexact2')
%     xlabel('x1')
%
%     figure(2)
%     plot(x1,error1,'k.-')
%     hold on
%     plot(x1,error2,'b.-')
%     title('Error')
%     legend('u1 error','u2 error')
%     xlabel('x1')
%
% %plots where the velocity is computed in a 2d domain
% else
%
%     figure(1)
%     plot(y1,y2,'k.-')
%     hold on
%     quiver(y1,y2,f1,f2,'r')
%     axis equal
%     quiver(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),u1m(1:skip:end,1:skip:end),u2m(1:skip:end,1:skip:end),'k')
%     xlim([x1min,x1max])
%     ylim([x2min,x2max])
%     title('Forces and Computed Velocity')
%
%     figure(2)
%     plot(y1,y2,'k.-')
%     hold on
%     quiver(y1,y2,f1,f2,'r')
%     axis equal
%     quiver(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),uexact1(1:skip:end,1:skip:end),uexact2(1:skip:end,1:skip:end),'k')
%     xlim([x1min,x1max])
%     ylim([x2min,x2max])
%     title('Forces and Exact Velocity')
%
%     figure(3)
%     plot(y1,y2,'k.-')
%     hold on
%     quiver(y1,y2,f1,f2,'r')
%     axis equal
%     quiver(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),error1(1:skip:end,1:skip:end),error2(1:skip:end,1:skip:end),5,'k')
%     xlim([x1min,x1max])
%     ylim([x2min,x2max])
%     title('Forces and Error in the Velocity')
%
%     figure(4)
%     pcolor(x1m,x2m,log10(errormag)+eps)
%     shading interp
%     colorbar
%     caxis( [-5 -3] )
%     hold on
%     plot(y1,y2,'k.-')
%     quiver(y1,y2,f1,f2,'r')
%     axis equal
%     xlim([x1min,x1max])
%     ylim([x2min,x2max])
%     title('Forces and Error Magnitude in the Velocity')
%
% end
%
%
%