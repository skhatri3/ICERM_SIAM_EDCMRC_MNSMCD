%main_permeable_cicle

%Based on Example 3 Section 3.1 of Cortez, Fluids, 2021
% Permeable circle with forces = curvature times normal vector
% Comppare to an exact solution
% Start with a circle of radius 1
% Compute the velocit from forces

% Developed for ICERM 
% July 30, 2024
% Author: Wanda Strychalski

close all
clear all
clc


video_flag = 0; % save the movie
if(video_flag)
    vidObj = VideoWriter('balloon_problem.mp4','MPEG-4');
    open(vidObj);
end


%% Parameters to set

%setting the viscosity
mu = 1;

% initial radius of the circle
r  = 1;

blob_num = 2;

%number of points on boundary where force is applied
N = 400; % match the paper

%discretization of cylinder boundary
ds = 2*pi*r/N; % initial
s = 0:ds:2*pi*r;
s = s(1:end-1)';

%regularization parameter
% ep = 2*ds;
ep = 0.1; % match the paper.
alpha = 1e-1; % was 1e-2
beta = 4/3*ep*alpha*ones(N,1);

%initial position cylinder on which forces are applied - cylinder is of radius r
y1 = r*cos(s/r);
y2 = r*sin(s/r);

y10 = r*cos(s/r); % initial configuration
y20 = r*sin(s/r);

% First derivative matrix without ds scaling
% periodic boundary conditions, centered differences
Dc = spdiags([-1 0 1],[-1:1],N,N);
Dc(1,N) = -1;
Dc(N,1) = 1;

% Second derivative matrix without ds scaling
% periodic boundary conditions, centered differences
DD2 = spdiags([1 -2 1],-1:1,N,N);
DD2(1,N) = 1;
DD2(N,1) = 1;

kappa = zeros(N,1); % initialize curvature

%resting radius of resulting circle
y1o = r*cos(s/r);
y2o = r*sin(s/r);

yo = [y1o, y2o];
dyo = Dc*yo;
magdyo = sqrt(dyo(:,1).^2 + dyo(:,2).^2);

%final time
tfinal = 1;
tfinal = r^2*mu/(2*alpha); % 5, the last time that the solution is valid

%number of timesteps (not specified in the paper)
tsteps = N*2;

%timestep
dtstep = tfinal/tsteps;
area = zeros(tsteps+1,1);
area(1) = polyarea([y1;y1(1)],[y2;y2(1)]);

time = 0:dtstep:tfinal;
time = time';
nt = length(time)-1;
Rcomputed = zeros(nt+1,1); % computed radius
Rexact = sqrt(r^2-2*alpha/mu*time); % exact solution

Rcomputed = zeros(tsteps,1);
R = 1/N*sum(sqrt(y1.^2 + y2.^2)); % formula from the paper
Rcomputed(1) = R;

dsArray = zeros(tsteps,1);
dsArray(1) = ds; % changes over time.

Rexact = sqrt(r^2-2*alpha/mu*time);
y = [y1, y2];

% In case we want to switch to Matlab's ODE solver later
% F = ode;
% F.Parameters = beta;
% F.InitialValue = y;
% F.ODEFcn = @(t,y,p) [y(2); p(1)*(1-y(1)^2)*y(2)-y(1)];
% F.Jacobian = @(t,y,p) [0 1; -2*p(1)*y(1)*y(2)-1  p(1)*(1-y(1)^2)];

%looping over time
for j=1:nt

    y = [y1, y2];

    %compute tangent
    D1 = 1/(2*ds).*Dc;
    dy = D1*y;
    magdy = sqrt(dy(:,1).^2 + dy(:,2).^2);
    tangent = dy./magdy;
    %compute normal
    normal = [tangent(:,2), -tangent(:,1)];

    D2 = 1/(ds^2)*DD2;
    dydy = D2*y;

    % compute curvature
    kappa = (dy(:,1).*dydy(:,2)-dy(:,2).*dydy(:,1))./(magdy.^3);

    beta = 4/3*ep*alpha*ones(N,1);

    f = (-dydy./magdy.^2);
    f = f*ds; % force and not force density
    b = beta.*f;

    %check_kappa = dot(g,normal,2);
    %diff = check_kappa-kappa;

    % quiver(y1,y2,normal(:,1),normal(:,2))
    % axis equal


    u_beta = RegStokeslets2D_permeable_gtovelocity(y,f,...
        y,ep,mu,blob_num, beta, normal);

    magu_beta = sqrt(u_beta(:,1).^2 + u_beta(:,2).^2);

    % u_exact = -alpha/mu*kappa.*normal;
    % magu_exact  = sqrt(u_exact(:,1).^2 + u_exact(:,2).^2);
    % diff = u_exact-u_beta;

    vdotn = dot(u_beta,normal,2); % for checking with the paper

    y1 = y1 + dtstep*(u_beta(:,1));
    y2 = y2 + dtstep*(u_beta(:,2));

    R = 1/N*sum(sqrt(y1.^2 + y2.^2));
    ds = 2*pi*R/N;
    dsArray(j+1) = ds;
    Rcomputed(j+1) = R;
    ep = 0.1*R; % dynamically update epsilon for accuracy

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

    %xlim([x1min,x1max])
    %ylim([x2min,x2max])
    axis equal
    time_str = sprintf('Time %.4f\n',time(j));
    title(time_str);
    hold off

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


    area(j) = polyarea([y1;y1(1)],[y2;y2(1)]);
end

if(video_flag)
    close(vidObj);
end


figure; % figure 2 in the paper (left)
semilogy(time,Rexact,time,Rcomputed,'--');
legend('Exact','Computed','Location','best');


figure; % figure 2 in the paper (right)
% last bit
ind = nt-150:nt;
plot(time(ind),Rexact(ind),time(ind),Rcomputed(ind),'--');
legend('Exact','Computed','Location','best');
grid on

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