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

% Create directory for saving figures
dformat = datetime; % Put date and time in directory name
dformat.Format = "ddMMMyyyy_'time_h'HH'm'mm's'ss";
directoryName = ['OutletVelocityError_StokesletsOnly_', char(dformat)];
mkdir(directoryName)
currPath = [pwd,'/',directoryName];
format long

%% Parameters to set 

% setting the viscosity
mu = 1; 
% Darcy number
Da = 0.3;

% number of points on boundary where velocity is set and force is computed 
Nvals = [40 80 120 160 200 240 280 320];
%[80 100 120 140 160 180 200 220 240 260 280 300 320];

%blob choice (1 is phi and 2 is psi)
blob=2;

% constant determining inflow velocity profile
G = 4;
% Scaled velocity profile at outlet
chi = 1;
 
% channel set-up
L = 3;                            % length of permeable region
c = 2;                            % length of solid region on half of channel
H = 1;                            % radius of channel
xmin = 0;
xmax = xmin + L + 2*c;
ymin = 0;
ymax = ymin + 2*H;
Ly = ymax - ymin;          % height of channel
Lx = xmax - xmin;          % length of channel
xm = (xmax + xmin)/2;  % midpoint of channel
perm_min = c;               % start of permeable part
perm_max = L + c;        % end of permeable part

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
epvals = logspace(-3,-0.4,20);

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
    %ds = (ymax-ymin)/N; % change later so that ds_x and ds_y are different
    ds_x = (xmax - xmin)/N;
    ds_y = (ymax - ymin)/( ceil((ymax - ymin)/ds_x));
    
    % usolutions=cell(1,length(N)); %for storing u solutions to compare
    % vsolutions=cell(1,length(N));

    uL2error = zeros(1,length(epvals));
    vL2error = zeros(1,length(epvals),1);
    uLinferror = zeros(1,length(epvals));
    vLinferror = zeros(1,length(epvals));

    for i = 1:length(epvals)
        ep1 = epvals(i);
        disp(['epsilon=', num2str(ep1)])

        %% Discretization
        %discretization of top/bottom walls
        stb = xmin+ds_x/2:ds_x:xmax-ds_x/2; stb = stb';
        %discretization of left/right wall (inlet/outlet)
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
        
        %coordinates of boundary source points
        y1 = [y1_top; y1_bot; y1_left]; %x-coordinates of all boundary points
        y2 = [y2_top; y2_bot; y2_left]; %y-coordinates of all boundary points
        
        %For use at very end:
        y1_right = xmax*ones(size(slr));   y2_right = slr; % right wall coordinates (y1_right,y2_right)
        %unit normals for right side
        normals_right=zeros(length(y1_right),2);
        normals_right(:,1)=1;
        far_idx = find( sqrt(ep1) < y2_right & y2_right < ymax-sqrt(ep1));
        y2_mid = y2_right(far_idx);
        y1_mid = xm*ones(size(y2_mid));

        %y1f=[y1; y1right];
        %y2f=[y2; y2right];

        % Define weights cooresponding to wall coordinates
        wt = [ds_x*ones(size(y1_top)); ds_x*ones(size(y1_bot)); ds_y*ones(size(y1_left))];

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
        u1_left=G*(y2_left/Ly.*(1-y2_left/Ly));
        u2_left=zeros(size(y1_left));
        
        u1 = [u1_top; u1_bot; u1_left];
        u2 = [u2_top; u2_bot; u2_left];

        %% MRS

        % Step 1: Use boundary velocities to compute the forces on the boundary
        f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep1,mu, blob,wt);
        f1 = f(:,1);
        f2 = f(:,2); 

        % Step 2: Use boundary forces to calculate velocities everywhere

        % Calculate velocities in middle
        %Umid=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[y1mid,y2mid],ep1,mu, blob, wt);
        %umid=Umid(:,1); vmid=Umid(:,2);

        % Calculate velocities at outlet
        Uout=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[y1_right,y2_right],ep1,mu, blob, wt);
        uout1=Uout(:,1); uout2=Uout(:,2);

        % Check flow rate at outlet 
        % inlet
        Rin = ds_y*sum(dot(normals_left, [u1_left u2_left]));

        %top
        % normals_top=zeros(length( y1_top),2);
        % normals_top(:,2)=1;
        %Rtop = ds_x*sum(dot(normals_top, [u1_top(1:length(y1_top)), u2_top(1:length(y1_top))]));

        idx_flux = find(y1_top >= c & y1_top <= xm);
        y1_flux = y1_top(idx_flux);
        normals_top_flux=zeros(length(y1_flux),2);
        normals_top_flux(:,2)=1;
        Rtop = ds_x*sum(dot(normals_top_flux, [u1_top(idx_flux), u2_top(idx_flux)]));
        
        % Change to evaluate at x = 2.5 instead of at outlet
        %x2_out= (ymin+ds/2:ds:ymax-ds/2)';
        %x1_out= xmax/2*ones(size(x2_out)); % xmax/2 for x = 2.5
        %uu = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],...
        %    [x1_out,x2_out],ep1,mu, blob);
        Rout = ds_y*sum(dot(normals_right, [uout1, uout2]));

        %% Calculate Current Errors
        % Using the difference between the inlet velocities and velocities
        % in the middle.  

        %diffu = Umid(:,1)-u1_side;
        %diffv = Umid(:,2)-u2_side;
        %diffu = Uout(:,1)-u1_mid;
        %diffv = Uout(:,2) - u2_mid;
        diffu = Uout(:,1) - u1_left;
        diffv = Uout(:,2) - u2_left;

        uL2error(i) = sqrt(ds_y*sum(diffu.^2))/norm(u1_left,2);
        vL2error(i) = sqrt(ds_x*sum(diffv.^2))/norm(u2_left,2);
        uLinferror(i) = max(max(abs(diffu)))/norm(u1_left,inf);
        vLinferror(i) = max(max(abs(diffv)))/norm(u2_left,inf);
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
    
 
end

%% Save plots

eppts = linspace(10^(-1.5),10^(-1),20);
yplot2 = 0.5 * eppts.^2;
%yplot2 = 4*eppts + 10^(-4);

figure(fig1); hold on;
%plot(eppts,yplot2,'--k','linewidth',2,'DisplayName','slope 2')
set(gca,'YScale','log');
set(gca,'xscale','log');
legend('Location','best')
box on
xlim([min(epvals), max(epvals)])
figsave_title = ['OutletVelocityError_L2_Permeable_StokesletsOnly','.fig'];
pngsave_title = ['OutletVelocityError_L2_Permeable_StokesletsOnly','.png'];
savefig([currPath,'/',figsave_title])
saveas(gcf, fullfile(currPath, pngsave_title))
hold off;

%fig_title = ['Horizontal Error (L2) - no permeability','.fig'];
%savefig([currPath,'/',fig_title])

figure(fig2)
hold on
%plot(eppts,yplot2,'--k','linewidth',2,'DisplayName','slope 2')
set(gca,'YScale','log');
set(gca,'xscale','log');
box on
xlim([min(epvals), max(epvals)])
legend('Location','best')
figsave_title = ['OutletVelocityError_Linf_Permeable_StokesletsOnly','.fig'];
pngsave_title = ['OutletVelocityError_Linf_Permeable_StokesletsOnly','.png'];
savefig([currPath,'/',figsave_title])
saveas(gcf, fullfile(currPath, pngsave_title))
hold off

% Save data for current run
clear fig1 fig2 % don't save figure handles
save(fullfile(currPath, 'RunData.mat'))

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
title('Best fit line for Epsilons w/minimal u-Error, Blob \psi')

legend('L2 Eps','L2 Line','Linf Eps','Linf Line', 'Location','best')

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
% legend('L2 Eps','L2 Line','Linf Eps','Linf Line', 'location', 'best')