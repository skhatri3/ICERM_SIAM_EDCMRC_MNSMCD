%% 2D Finite Channel, Nonpermeable, flow pics and rates

% Channel set-up based off of Example 3 of Cortez, Fluids 2021
% Channel has prescribed parabolic inflow on left, solid walls on top and
% bottom, and an open outlet on the right.

% Nonpermeable Channel: Using this set-up to better understand the
% convergence behavior of the Stokeslets in the 2D, permeable, finite
% channel problem.

% Looking at channel without permeability, examining near and far field
% errors. We get first-order for far-field and half-order for near. This
% is expected since we have corners in our domain.

clear all 

%currDate = datestr(datetime);

%mkdir(currDate)
%currPath = [pwd,'/',currDate];
format long

%% Parameters to set 

% setting the viscosity
mu = 1; 

% number of points on boundary where velocity is set and force is computed 
Nvals = [40 80 120 160 200 240 280 320];
%[80 100 120 140 160 180 200 220 240 260 280 300 320];


% % C1 value range
% C1vals = linspace(0.1,4,20);

% % Keeping ep2 constant (looking for C1 in ep1 = C1ds)
% %ep2 = 0.0002;
% ep2 = 0.05/sqrt(5);

%blob choice (1 is phi and 2 is psi)
blob=2;

%constant determining inflow velocity profile
a=4;
 
%channel set-up
Lx = 5;
Ly = 2;
xmin = 0;
xmax = xmin + Lx;
ymin = 0;
ymax = ymin + Ly;
perm_min=5/3;  %start of permeable part
perm_max=10/3;  %end of permeable part


% Initialize Plots
fig1 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Horizontal Velocity Error (L2), non-permeable, using blob $\phi$',...
    'FontSize', 18,'Interpreter','latex')
hold off;


fig2 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Horizontal Velocity Error (Linf), non-permeable, using blob $\phi$',...
    'FontSize', 18,'Interpreter','latex') 
hold off;

 fig3 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Vertical Velocity Error (L2), non-permeable, using blob $\psi$',...
    'FontSize', 18,'Interpreter','latex')
hold off;


fig4 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Vertical Velocity Error (Linf), non-permeable, using blob $\psi$',...
    'FontSize', 18,'Interpreter','latex') 
hold off;

% epvals=0.00001*2.^(0:16);
% epvals=[epvals (0.005:0.001:0.2)];
% epvals=sort(epvals);

epvals = logspace(-3,-0.5,20);
% Set up error matrix and initialize figure
uErrorMatrixL2 = zeros(length(Nvals),length(epvals));
vErrorMatrixL2 = uErrorMatrixL2;
uErrorMatrixLinf = uErrorMatrixL2;
vErrorMatrixLinf = uErrorMatrixL2;

%%
for k = 1:length(Nvals)
    N = Nvals(k)
    %discretization of channel
    ds = (ymax-ymin)/N;
    
    % usolutions=cell(1,length(N)); %for storing u solutions to compare
    % vsolutions=cell(1,length(N));

    uL2error = zeros(1,length(epvals));
    vL2error = zeros(1,length(epvals),1);
    uLinferror = zeros(1,length(epvals));
    vLinferror = zeros(1,length(epvals));

    for i = 1:length(epvals)
        ep1 = epvals(i)

        %% Discretization

        %discretization of top and bottom:
        s = ds/2:ds:xmax-ds/2;
        s = s';
        %top wall coordinates (x,y) = (y1_top,y2_top)
        y1_top=s;
        y2_top = ymax*ones(size(y1_top));
        %unit normals for top: 
        normals_top=zeros(length(y1_top),2);
        normals_top(:,2)=1;
        %bottom wall coordinates (x,y) = (y1_bot,y2_bot)
        y1_bot = s;
        y2_bot = ymin*ones(size(y1_bot));
        %unit normals for bottom
        normals_bot=zeros(length(y1_bot),2); 
        normals_bot(:,2)=-1;
        %left-hand wall coordinates (x,y) = (y1_side,y2_side)
        y2_side = (ymin+ds/2:ds:ymax-ds/2)';
        y1_side = xmin*ones(size(y2_side));
        %unit normals for side
        normals_side=zeros(length(y1_side),2);
        normals_side(:,1)=-1;

        %coordinates of boundary points
        y1 = [y1_top;y1_bot; y1_side];
        y2 = [y2_top;y2_bot; y2_side];
        %For use at very end:
        y2right = (ymin+ds/2:ds:ymax-ds/2)';
        y1right = xmax*ones(size(y2_side));

        %y2mid = y2right;
        %y1mid = xmax/2*ones(size(y2mid));

        far_idx = find( sqrt(ep1) < y2right & y2right < Ly-sqrt(ep1));
        y2mid = y2right(far_idx);
        y1mid = xmax/2*ones(size(y2mid));

        y1f=[y1; y1right];
        y2f=[y2; y2right];

        % Define weights cooresponding to wall coordinates
        wt = ds*ones(size(y1));

        % Setting velocities on boundary
        % top - no flow
        u1_top=zeros(size(y1_top));
        u2_top=zeros(size(y1_top));
        % bottom - no flow
        u1_bot=zeros(size(y1_bot));
        u2_bot=zeros(size(y1_bot));
        % inlet: poiseuille flow
        u1_side=a*(y2_side/Ly.*(1-y2_side/Ly));
        u2_side=zeros(size(y1_side));
        u1_mid = a*(y2mid/Ly.*(1-y2mid/Ly));
        u2_mid = zeros(size(y1mid));

        u1 = [u1_top; u1_bot; u1_side];
        u2 = [u2_top; u2_bot; u2_side];

        %% MRS

        % Step 1: Use boundary velocities to compute the forces on the boundary
        f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep1,mu, blob,wt);
        f1 = f(:,1);
        f2 = f(:,2); 

        % Step 2: Use boundary forces to calculate velocities everywhere

        % Calculate velocities in middle
        Umid=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[y1mid,y2mid],ep1,mu, blob, wt);
        umid=Umid(:,1); vmid=Umid(:,2);


        % Check flow rates 

        % inlet
        % Rin=ds*sum(dot(normals_side, [u1_side u2_side]));

        %top
        % normals_top=zeros(length( y1_top),2);
        % normals_top(:,2)=1;
        % Rtop=ds*sum(dot(normals_top, [ufull(1:length(y1_top)), vfull(1:length(y1_top))]));

        % Change to evaluate at x = 2.5 instead of at outlet
        %x2_out= (ymin+ds/2:ds:ymax-ds/2)';
        %x1_out= xmax/2*ones(size(x2_out)); % xmax/2 for x = 2.5
        %uu = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],...
        %    [x1_out,x2_out],ep1,mu, blob);
        %Rout=-ds*sum(dot(normals_side, uu));

        %% Calculate Current Errors
        % Using the difference between the inlet velocities and velocities
        % in the middle.  

        %diffu = Umid(:,1)-u1_side;
        %diffv = Umid(:,2)-u2_side;
        diffu = Umid(:,1)-u1_mid;
        diffv = Umid(:,2) - u2_mid;

        uL2error(i) = sqrt(ds*sum(diffu.^2))/norm(u1_mid,2);
        vL2error(i) = sqrt(ds*sum(diffv.^2))/norm(u2_mid,2);
        uLinferror(i) = max(max(abs(diffu)))/norm(u1_mid,inf);
        vLinferror(i) = max(max(abs(diffv)))/norm(u2_mid,inf);
        %error(i) = norm()
        %error(i)=Rin+Rtop+Rout;

    end

    %%
    % Save Current Error results and update figure
    uErrorMatrixL2(k,:) = uL2error; % save in row k
    uErrorMatrixLinf(k,:) = uLinferror;
    vErrorMatrixL2(k,:) = vL2error; % save in row k
    vErrorMatrixLinf(k,:) = vLinferror;


    name = ['N = ' num2str(N)];
    
    figure(fig1); hold on;
    loglog(epvals, uL2error,'o-', 'LineWidth', 2.5, 'MarkerSize', ...
       10,'DisplayName',name)
    hold off;

    figure(fig2); hold on;
    loglog(epvals,uLinferror,'x--','Linewidth',2.5, 'MarkerSize',10,'DisplayName',name)
    hold off;
    
     figure(fig3); hold on;
     loglog(epvals, vL2error,'o-', 'LineWidth', 2.5, 'MarkerSize', ...
        10,'DisplayName',name)
     hold off;
    % 
    % figure(fig4); hold on;
    % loglog(epvals,vLinferror,'x--','Linewidth',2.5, 'MarkerSize',10,'DisplayName',name)
    % hold off;

 
end

%% Save plots

eppts = linspace(10^(-1.5),10^(-1),20);
yplot2 = 0.5 * eppts.^2;
%yplot2 = 4*eppts + 10^(-4);

figure(fig1); hold on;
%plot(eppts,yplot2,'--k','linewidth',2,'DisplayName','slope 2')
set(gca,'YScale','log');
set(gca,'xscale','log');
legend
hold off;

%fig_title = ['Horizontal Error (L2) - no permeability','.fig'];
%savefig([currPath,'/',fig_title])

figure(fig2)
hold on
%plot(eppts,yplot2,'--k','linewidth',2,'DisplayName','slope 2')
set(gca,'YScale','log');
set(gca,'xscale','log');
legend


%fig_title = ['Horizontal Error (Linf)- no permeability','.fig'];
%savefig([currPath,'/',fig_title])

%% Best fit lines

% u-Error
[minL2,indexL2] = min(uErrorMatrixL2(2:end,:),[],2); % minimize along rows, since rows correspond to fixed ds
[minLinf,indexLinf] = min(uErrorMatrixLinf(2:end,:),[],2);
figure
hold on
ds = 1./Nvals(2:end);
plot(ds,epvals(indexL2),'ro','linewidth',2,'markersize',10)
p2 = polyfit(ds,epvals(indexL2),1);
q2 = polyfit(epvals(indexL2),uErrorMatrixL2(indexL2),1);
l2 = polyval(p2,ds);
plot(ds,l2,'--r','linewidth',2)

plot(ds,epvals(indexLinf),'bx','linewidth',2,'markersize',10)
pinf = polyfit(ds,epvals(indexLinf),1);
linf = polyval(pinf,ds);
plot(ds,linf,'-.b','linewidth',2)

xlabel('ds=1/N'), ylabel('Epsilon')
title('Best fit line for Epsilons w/minimal u-Error, Blob \phi')

legend('L2 Eps','L2 Line','Linf Eps','Linf Line')

%fig_title = ['Best_Fit_Line','.fig'];
%savefig([currPath,'/',fig_title])

%% v-Error
% [vminL2,vindexL2] = min(vErrorMatrixL2,[],2); % minimize along rows, since rows correspond to fixed ds
% [vminLinf,vindexLinf] = min(vErrorMatrixLinf,[],2);
% figure
% hold on
% plot(ds,epvals(vindexL2),'ro','linewidth',2,'markersize',10)
% vp2 = polyfit(ds,epvals(vindexL2),1);
% vl2 = polyval(vp2,ds);
% plot(ds,vl2,'--r','linewidth',2)
% 
% plot(ds,epvals(vindexLinf),'bx','linewidth',2,'markersize',10)
% vpinf = polyfit(ds,epvals(vindexLinf),1);
% vlinf = polyval(vpinf,ds);
% plot(ds,vlinf,'-.b','linewidth',2)
% 
% xlabel('ds=1/N'), ylabel('Epsilon')
% title('Best fit line for Epsilons w/minimal v-Error')
% 
% legend('L2 Eps','L2 Line','Linf Eps','Linf Line')