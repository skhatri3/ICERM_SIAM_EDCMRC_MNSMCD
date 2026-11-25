%% 2D Channel using Pozrikidis geometry with permeable portion of length L, flow pics and rates

% Channel set-up based off of Example 3 of Cortez, Fluids 2021
% Channel has prescribed parabolic inflow on left, solid walls on top and
% bottom, and an open outlet on the right.

% Channel with transmural pressure boundary condition: Using this set-up to better understand the
% convergence behavior of the Stokeslets in the 2D, permeable, channel problem.

% Looking at channel with stokeslets only, finding
% relationship between ds and epsilon using L2 and Linf errors.
% Later will apply ideas using MRS+SD

% Modified by Kurianski for Pozrikidis channel Nov. 2026

clear; close all;

addpath('../stokeslets_codes'); 
%currDate = datestr(datetime);

%mkdir(currDate)
%currPath = [pwd,'/',currDate];
format long

%% Parameters to set 

% setting the viscosity
mu = 1;
% Darcy number
Da = 0.3;

% number of points on boundary where velocity is set and force is computed 
Nvals = [40 60 80 100 120 140 160 180 200 220 240 260 280 300 320];

% blob choice (1 is phi and 2 is psi)
blob = 2;
% constant determining inflow velocity profile
G = 4;
% Scaled velocity profile at outlet
chi = 1;
 
% channel set-up
L = 3;
c = 2;
H = 1;
xmin = 0;
xmax = xmin + L + 2*c;
ymin = 0;
ymax = ymin + 2*H;
perm_min = c;  %start of permeable part
perm_max = L + c;  %end of permeable part
Ly = ymax - ymin;
Lx = xmax - xmin;
xm = (xmax + xmin)/2; % midpoint of channel

% Use Newton's method to compute eta eigenvalue
eta = sqrt(2*Da);
fun = @(lam) (Da-H/2)*lam*cot(lam*H)*cos(lam*H) + cos(lam*H)/2 - H*lam*sin(lam*H)/2;
eta = fzero(fun, eta);

% Initialize Plots
fig1 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Horizontal Velocity Error (L2), Stokeslets Only, using blob $\psi$',...
    'FontSize', 18,'Interpreter','latex')
hold off;


fig2 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Horizontal Velocity Error (Linf), Stokeslets Only, using blob $\psi$',...
    'FontSize', 18,'Interpreter','latex') 
hold off;

% Epsilon values
epvals = logspace(-3,-0.5,20);

% Set up error matrix and initialize figure
uErrorMatrixL2 = zeros(length(Nvals),length(epvals));
vErrorMatrixL2 = uErrorMatrixL2;
uErrorMatrixLinf = uErrorMatrixL2;
vErrorMatrixLinf = uErrorMatrixL2;

%%
for k = 1:length(Nvals)
    N = Nvals(k);
    disp(['N=', num2str(N)]);
    %discretization of channel
    ds_x = (xmax - xmin)/N;
    ds_y = (ymax - ymin)/( ceil((ymax - ymin)/ds_x));
    
    % usolutions=cell(1,length(N)); %for storing u solutions to compare
    % vsolutions=cell(1,length(N));

    uL2error = zeros(1,length(epvals));
    vL2error = zeros(1,length(epvals));
    uLinferror = zeros(1,length(epvals));
    vLinferror = zeros(1,length(epvals));
    flowRateError = zeros(1,length(epvals));

    for i = 1:length(epvals)
        ep1 = epvals(i);
        disp(['epsilon=', num2str(ep1)])

        %% Discretization
        %discretization step of top/bottom walls
        stb = xmin+ds_x/2:ds_x:xmax-ds_x/2; stb = stb';
        %discretization step of left/right wall (inlet/outlet)
        slr = ymin+ds_y/2:ds_y:ymax-ds_y/2; slr = slr';
        % Define coordinates on each wall (x-coord,y-coord) (source points)
        y1_top = stb;   y2_top = ymax*ones(size(stb)); % top wall coordinates (y1_top,y2_top)
        %unit normals for top:
        normals_top=zeros(length(y1_top),2);
        normals_top(:,2)=1;
        y1_bot = stb;   y2_bot = ymin*ones(size(stb)); % bottom wall coordinates (y1_bot,y2_bot)
        %unit normals for bottom
        normals_bot=zeros(length(y1_bot),2);
        normals_bot(:,2)=-1;
        y1_left = xmin*ones(size(slr));   y2_left = slr; % left wall coordinates (y1_left,y2_left)
        %unit normals for left side
        normals_left=zeros(length(y1_left),2);
        normals_left(:,1)=-1;
        y1_right = xmax*ones(size(slr));   y2_right = slr; % right wall coordinates (y1_right,y2_right)
        %unit normals for right side
        normals_right=zeros(length(y1_right),2);
        normals_right(:,1)=1;
        y1 = [y1_top; y1_bot; y1_left; y1_right]; %x-coordinates of all boundary points
        y2 = [y2_top; y2_bot; y2_left; y2_right]; %y-coordinates of all boundary points
        y = [y1, y2];

        far_idx = find( ep1 < y2_right & y2_right < Ly-ep1);
        y2_mid = y2_right(far_idx);
        y1_mid = xmax/2*ones(size(y2_mid));

        % Define weights cooresponding to wall coordinates using midpoint rule
        wt = [ds_x*ones(size(y1_top)); ds_x*ones(size(y1_bot)); ds_y*ones(size(y1_left)); ds_y*ones(size(y1_right))];

        % Setting velocities on boundary
        % top velocity (no flow at solid walls; v = Da*p_T in permeable)
        u1_top=zeros(size(y1_top));
        u2_top=zeros(size(y1_top));
        % bottom velocity (no flow at solid walls; v = -Da*p_T in permeable)
        u1_bot=zeros(size(y1_bot));
        u2_bot=zeros(size(y1_bot));
        % Analytical velocity in permeable region
        idx_perm = find(perm_min <= y1_top & y1_top <= perm_max);
        u2_top(idx_perm) = -Da*G*sinh(eta*(y1_top(idx_perm)-xm))/(eta*cosh(eta*L/2));
        u2_bot(idx_perm) = Da*G*sinh(eta*(y1_top(idx_perm)-xm))/(eta*cosh(eta*L/2));
        % inlet: poiseuille flow
        u1_left=G*(y2_left/(2*H)).*(1-y2_left/(2*H));
        u2_left=zeros(size(y1_left));
        % outlet: scaled poiseuille flow
        u1_right = chi*G*(y2_left/(2*H)).*(1-y2_left/(2*H));
        u2_right = zeros(size(y1_left));
        
        u1 = [u1_top; u1_bot; u1_left; u1_right];
        u2 = [u2_top; u2_bot; u2_left; u2_right];

        %% MRS
        % Step 1: Use boundary velocities to compute the forces on the boundary
        f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep1,mu, blob,wt);
        f1 = f(:,1);
        f2 = f(:,2); 

        % Step 2: Use boundary forces to calculate velocities everywhere
        % Calculate velocities in middle
        Umid=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[y1_mid,y2_mid],ep1,mu, blob, wt);
        umid=Umid(:,1); vmid=Umid(:,2);

        % Check flow rates 
        % inlet flow rate
        Rin = ds_x*sum(dot(normals_left, [u1_left, u2_left]));
        % top flow rate
        Rtop = ds_y*sum(dot(normals_top, [ufull(1:length(y1_top)), vfull(1:length(y1_top))]));
        % outlet flow rate
        Rout = ds_x*sum(dot(normals_right, [u1_right, u2_right]));

        %% Calculate Current Errors
        % Using the difference between the inlet velocities and velocities
        % in the middle.  
        diffu = Umid(:,1)-u1_mid;
        diffv = Umid(:,2) - u2_mid;

        %uL2error(i) = sqrt(ds_y*sum(diffu.^2))/norm(u1_mid,2);
        %vL2error(i) = sqrt(ds_y*sum(diffv.^2))/norm(u2_mid,2);
        %uLinferror(i) = max(max(abs(diffu)))/norm(u1_mid,inf);
        %vLinferror(i) = max(max(abs(diffv)))/norm(u2_mid,inf);
        flowRateError(i) = Rin+Rtop+Rout;

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
    
 
end

%% Save plots

eppts = linspace(10^(-1.5),10^(-1),20);
yplot2 = 0.5 * eppts.^2;
%yplot2 = 4*eppts + 10^(-4);

figure(fig1); hold on;
plot(eppts,yplot2,'--k','linewidth',2,'DisplayName','slope 2')
set(gca,'YScale','log');
set(gca,'xscale','log');
legend
hold off;

%fig_title = ['Horizontal Error (L2) - no permeability','.fig'];
%savefig([currPath,'/',fig_title])

figure(fig2)
set(gca,'YScale','log');
set(gca,'xscale','log');
legend


%fig_title = ['Horizontal Error (Linf)- no permeability','.fig'];
%savefig([currPath,'/',fig_title])

%% Best fit lines

% u-Error
[minL2,indexL2] = min(uErrorMatrixL2,[],2); % minimize along rows, since rows correspond to fixed ds
[minLinf,indexLinf] = min(uErrorMatrixLinf,[],2);
figure
hold on
ds_y = 1./Nvals;
plot(ds_y,epvals(indexL2),'ro','linewidth',2,'markersize',10)
p2 = polyfit(ds_y,epvals(indexL2),1);
q2 = polyfit(epvals(indexL2),uErrorMatrixL2(indexL2),1);
l2 = polyval(p2,ds_y);
plot(ds_y,l2,'--r','linewidth',2)

plot(ds_y,epvals(indexLinf),'bx','linewidth',2,'markersize',10)
pinf = polyfit(ds_y,epvals(indexLinf),1);
linf = polyval(pinf,ds_y);
plot(ds_y,linf,'-.b','linewidth',2)

xlabel('ds=1/N'), ylabel('Epsilon')
title('Best fit line for Epsilons w/minimal u-Error, Blob \psi')

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