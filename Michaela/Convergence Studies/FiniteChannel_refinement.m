%Example 3 of Cortez, Fluids 2021
%Channel with inflow and part of membrane permeable

%Developed by Ricardo Cortez, Brittany Leathers, and Michaela Kubacki
%July 2024

% Using 2 different epsilons, calculating convergence rates over a grid
% that includes a boundary buffer of 3*ep.


clear all 
close all

format long

%% Parameters to set 

%setting the viscosity
mu = 1; 

%number of points on boundary where velocity is set and force is computed 
%Nvals = [15, 30, 60, 120, 240, 480];
Nvals = [10,20,40,80,160,320,640];%,640];


%blob choice
blob=2;

%constant determining inflow velocity profile
a=4;

 
%channel
Lx = 5;
Ly = 1;
xmin = 0;
xmax = xmin + Lx;
ymin = 0;
ymax = ymin + Ly;
perm_min=5/3;  %start of permeable part
perm_max=10/3;  %end of permeable part


usolutions=cell(1,length(Nvals)); %for storing u solutions to compare
vsolutions=cell(1,length(Nvals));
error = zeros(1,length(Nvals));
errorOut = zeros(1,length(Nvals));
errorTop = zeros(1,length(Nvals));


%%
for i=1:length(Nvals)
    N=Nvals(i)

    %discretization of channel
    ds = (ymax-ymin)/N;

    %%
    %regularization parameters
    % C-values: Using values determined numerically 
    %C1 = 3.5;%1.665;%3.5;% 3.3;
    C2 = 0.05/sqrt(5)*(160)^(2/5);%0.05/sqrt(5)*(160)^(1/3);%0.34;%0.06;%0.012747
    C1 = 0.05/sqrt(5)*(160)^(1/2);
    ep1 = C1*ds.^(1/2);
    ep2 = C2*ds.^(1/5);
    %ep1 = ep2;

% Previously used values when only used 1 epsilon:
%ep = 0.05/sqrt(5)*160^(1/3)*ds^(1/3);
%ep = 0.05/sqrt(5)*160^(1/4)*ds^(1/4); 

%permeability coefficient (assuming constant for the permeable region) 
% for blob 1
%b = 0.011843080738230*ep; % for b = 0.0001 (beta 1)
%b = 0.004247164369662*ep; % for b = 0.00018
%b = 0.0282306994128518*ep; % for b = 0.00024 flow


%permeability coefficient (assuming constant for the permeable region) 
% for blob 2
% b=0.010733*ep; % for b = 0.00024 (beta 3)
% b = 2.028e-07 - 0.0001483*ep1 + 0.008181*ep2;
b = 0.008136682703766*ep2 - 0.000001993238909; % for b = 0.00018 (beta 2)
% b =0.0045233552698796*ep2 - 0.000001179437583; % for b = 0.0001 (beta 1)

%discretization of top and bottom (using midpoints):
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
y1f=[y1; y1right];
y2f=[y2; y2right];

% Define weights cooresponding to wall coordinates
wt = ds*ones(size(y1));

%indices of permeable region
I=find(y1_top<2/3*Lx & y1_top>1/3*Lx);
I2=find(y1_top>=2/3*Lx | y1_top<=1/3*Lx);
%unit normals
normals=[normals_top; normals_bot; normals_side];


%beta vector for function;
beta= zeros(length(y1),1); 
%beta(I)=b;



%velocity on boundary
%top
u1_top=zeros(size(y1_top));
u2_top=zeros(size(y1_top));
%bottom
u1_bot=zeros(size(y1_bot));
u2_bot=zeros(size(y1_bot));
%side: poiseulle flow
u1_side=a*(y2_side/Ly.*(1-y2_side/Ly));
u2_side=zeros(size(y1_side));

u1 = [u1_top; u1_bot; u1_side];
u2 = [u2_top; u2_bot; u2_side];


%
% Step 1: Compute g on boundary (using St and SD)
g=RegStokeslets2D_velocityto_gforce_permeable_diff_ep([y1,y2],[y1,y2],...
    [u1,u2], ep1, ep2, mu, blob, I, beta, normals,wt);

% CHANGE 1:
%g = g/ds;

% Step 2: Find (missing) velocity in permeable region (due to SD only)
y1b=y1(I);
y2b=y2(I);
[u_beta]=RegStokeslets2D_gtovelocity([y1,y2],g, [y1b,y2b],...
    ep2,mu, blob, beta, normals,wt);

% CHANGE 2:
%u_beta = u_beta*ds;

%Put in the new velocities for the permeable part
%permeable part : temporary velocity for finding g
u1(I)=u_beta(:,1);
u2(I)=u_beta(:,2);

%Step 3: Compute Forces everywhere on boundary (St) 
f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep1,mu, blob,wt);
% CHANGE 3:
%f = f/ds;
f1 = f(:,1);
f2 = f(:,2);



% Step 4: Compute flow in the channel (St)
% Calculate velocities on right hand side too
Ufull=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[y1f,y2f],ep1,mu, blob, wt);
% CHANGE 4:
%Ufull = Ufull*ds;

ufull=Ufull(:,1); vfull=Ufull(:,2);




%% Calculate on Grid
% Need a buffer from the boundary for optimal convergence. Using 3ep.

dist = 0.1;
%dist = 3*(C2*(1/Nvals(1))^(1/3)); % buffer from boundary using largest ep2
%dist = 3*C1*(1/Nvals(1))% buffer from boundary using largest ep1
dxg = 1/200;
Nx = round((xmax-2*dist)/dxg);
Ny = round((ymax - 2*dist)/dxg);
xg = linspace(dist,xmax-dist,Nx);
yg = linspace(dist,ymax-dist,Ny);
[xgg,ygg] = meshgrid(xg,yg);
xg = xgg(:);
yg = ygg(:);
Ugrid=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[xg,yg],ep1,mu,blob, wt);
% CHANGE 5
%Ugrid = Ugrid*ds;
ug=reshape(Ugrid(:,1), Nx, Ny); 
vg=reshape(Ugrid(:,2), Nx, Ny); 
speed=sqrt(ug.^2+vg.^2);
usolutions{i} = ug;
vsolutions{i} = vg;
xgo = xg';
ygo = yg';

%% Check flow rates 

% inlet
Rin=ds*sum(dot(normals_side, [u1_side u2_side]));

%top
normals_top=zeros(length( y1_top),2);
normals_top(:,2)=1;
Rtop=ds*sum(dot(normals_top, [ufull(1:length(y1_top)), vfull(1:length(y1_top))]));

%outlet
x2_out= (ymin+ds/2:ds:ymax-ds/2)';
x1_out= xmax*ones(size(x2_out));
uu = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],...
    [x1_out,x2_out],ep1,mu, blob,wt);
% CHANGE 6
%uu = uu*ds;
Rout=-ds*sum(dot(normals_side, uu));

error(i)=Rin+Rtop+Rout;
errorOut(i)=Rout;
errorTop(i) = Rtop;

end
%% Net Flow Plot

colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

% for plotting 1/N and 1/N^2 in error plots
forploty0=0.1;
forplotdx0=20;
dxforplot=1./Nvals;
yforplot=forploty0/forplotdx0^1*dxforplot.^1;
forploty02=20;
forplotdx02=10;
dxforplot2=1./Nvals;
yforplot2=forploty02/forplotdx02^2*dxforplot2.^2;
   
% Net Flow Error Plot
figure;
loglog(Nvals, abs(error),'o-', 'LineWidth', 2.5, 'MarkerSize', ...
       10, 'Color', colorlb)
hold on
loglog(Nvals, yforplot2, 'LineWidth', 2,'Color', colorg)
text(200,1.6/1000000,'$1/N^2$', 'Color', colorg, 'FontSize', 14, ...
    'Interpreter', 'latex');
hold off
axis([min(Nvals) max(Nvals) 1e-10 1]);
set(gca, 'FontSize', 17);
xlabel('$N$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Net Flow Rate, $\beta = 0.00018$, $\epsilon_1 = Cds^{1/2}$, $\epsilon_2 = C(ds)^{2/5}$, blob $\psi$',...
    'FontSize', 18,'Interpreter','latex') 
pError = [0 abs(error(1:end-1))./abs(error(2:end))];

table(Nvals',abs(error)',pError','VariableNames',{'N','Net Flow Error','Order-Est'})

% Rout Error Plot
% figure;
% %loglog(Nvals, abs(errorOut),'o-', 'LineWidth', 2.5, 'MarkerSize', ...
% %       10, 'Color', colorlb)
% hold on
% loglog(Nvals, yforplot2, 'LineWidth', 2,'Color', colorg)
% text(200,1.6/1000000,'$1/N^2$', 'Color', colorg, 'FontSize', 14, ...
%     'Interpreter', 'latex');
% hold off
% axis([min(Nvals) max(Nvals) 1e-10 1]);
% set(gca, 'FontSize', 17);
% xlabel('$N$', 'FontSize', 18,'Interpreter','latex')
% ylabel('Error', 'FontSize', 17)
% title('Outflow Error, $\beta = 0.00018$, $\epsilon_1 = C ds^{1/2}$, $\epsilon_2 = C(ds)^{2/5}$, blob $\psi$',...
%     'FontSize', 18,'Interpreter','latex') 
% pOutError = [0 abs(errorOut(1:end-1))./abs(errorOut(2:end))];

%table(Nvals',abs(errorOut)',pOutError','VariableNames',{'N','Out Flow Error','Order-Est'})

% Rtop Error Plot
% figure;
% %loglog(Nvals, abs(errorTop),'o-', 'LineWidth', 2.5, 'MarkerSize', ...
% %       10, 'Color', colorlb)
% hold on
% %loglog(Nvals, yforplot2, 'LineWidth', 2,'Color', colorg)
% text(200,1.6/1000000,'$1/N^2$', 'Color', colorg, 'FontSize', 14, ...
%     'Interpreter', 'latex');
% hold off
% axis([min(Nvals) max(Nvals) 1e-10 1]);
% set(gca, 'FontSize', 17);
% xlabel('$N$', 'FontSize', 18,'Interpreter','latex')
% ylabel('Error', 'FontSize', 17)
% title('Topflow Error, $\beta = 0.00018$, $\epsilon_1 = C ds^{1/2}$, $\epsilon_2 = C(ds)^{2/5}$, blob $\psi$',...
%     'FontSize', 18,'Interpreter','latex') 
% pTopError = [0 abs(errorTop(1:end-1))./abs(errorTop(2:end))];

%table(Nvals',abs(errorTop)',pOutError','VariableNames',{'N','Top Flow Error','Order-Est'})




%% Calculate successive difference errors
d1u=zeros(1, length(Nvals));
d2u=zeros(1, length(Nvals));
dmaxu=zeros(1, length(Nvals));

d1v=zeros(1, length(Nvals));
d2v=zeros(1, length(Nvals));
dmaxv=zeros(1, length(Nvals));

for i=1:length(Nvals)-1
    % calculate differences from fine to coarse 
    ufine=usolutions{i+1}; 
    ucoarse=usolutions{i};
    diffu=(ucoarse-ufine);
    vfine=vsolutions{i+1}; 
    vcoarse=vsolutions{i};
    diffv=(vcoarse-vfine);    
    % compute difference norms L1, L2, Linf
    dx=dxg;
    d1u(1, i)= dx^2*sum(sum(abs(diffu)));
    d2u(1, i)= sqrt(dx^2*sum(sum(diffu(:,:).^2)));
    dmaxu(1, i)= max(max(abs(diffu)));
    d1v(1, i)= dx^2*sum(sum(abs(diffv)));
    d2v(1, i)= sqrt(dx^2*sum(sum(diffv(:,:).^2)));
    dmaxv(1, i)= max(max(abs(diffv)));

end
%% Tables of errors with order estimates
p1u = [0 d1u(1,1:end-2)./d1u(1,2:end-1)];
p1v = [0 d1v(1,1:end-2)./d1v(1,2:end-1)];
p2u = [0 d2u(1,1:end-2)./d2u(1,2:end-1)];
p2v = [0 d2v(1,1:end-2)./d2v(1,2:end-1)];
pmaxu = [0 dmaxu(1,1:end-2)./dmaxu(1,2:end-1)];
pmaxv = [0 dmaxv(1,1:end-2)./dmaxv(1,2:end-1)];
table(Nvals(2:end)',d1u(1,1:end-1)',p1u',d1v(1,1:end-1)',p1v','VariableNames',{'N','L1-u','u_orderEst','L1-v','v_orderEst'})
table(Nvals(2:end)',d2u(1,1:end-1)',p2u',d2v(1,1:end-1)',p2v','VariableNames',{'N','L2-u','u_orderEst','L2-v','v_orderEst'})
table(Nvals(2:end)',dmaxu(1,1:end-1)',pmaxu',dmaxv(1,1:end-1)',pmaxv','VariableNames',{'N','Linf-u','u_orderEst','Linf-v','v-orderEst'})

%% Allnorms refinement plots
   
colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

forploty0=100;
forplotdx0=10;
dxforplot=1./Nvals;
yforplot=forploty0/forplotdx0^1*dxforplot.^1;
forploty02=100;
forplotdx02=10;
dxforplot2=1./Nvals;
yforplot2=forploty02/forplotdx02^2*dxforplot2.^2;

% Horizontal Velocity

figure;
% Inf-Norm
loglog(Nvals, dmaxu,'o-', 'LineWidth', 3, 'MarkerSize', ...
       10, 'Color', colorlb)
hold on
% 2-Norm
loglog(Nvals, d2u, 's-', 'LineWidth', 3,'MarkerSize', 10, 'Color', colordb);
% 1-Norm
loglog(Nvals, d1u, 'd-', 'LineWidth', 3,'MarkerSize', 10,'Color', colorp);
% Compare to 1/N^2
loglog(Nvals, yforplot2, 'LineWidth', 2,'Color', colorg)
text(100,1/10000,'1/N^2', 'Color', colorg, 'FontSize', 16);
% Compare to 1/N
loglog(Nvals, yforplot,'--','LineWidth', 2,'Color', colorg)
text(100,1/10000,'1/N', 'Color', colorg, 'FontSize', 16);
hold off
legend('L^{\infty}', 'L^2', 'L^1');
axis([min(Nvals) Nvals(end-1) 10^(-6) 10]);
set(gca, 'FontSize', 12);
xlabel('N', 'FontSize', 18)
ylabel('Error', 'FontSize', 18)
title('u-Error,  $\beta = 0.00018$, $\epsilon_1 = C ds^{1/2}$, $\epsilon_2 = C(ds)^{2/5}$, using blob $\psi$',...
    'FontSize', 18,'Interpreter','latex')
   

% Vertical Velocity
figure;
% Inf-Norm
loglog(Nvals, dmaxv,'o-', 'LineWidth', 3, 'MarkerSize', ...
       10, 'Color', colorlb)
hold on
% 2-Norm
loglog(Nvals, d2v, 's-', 'LineWidth', 3,'MarkerSize', 10, 'Color', colordb);
% 1-Norm
loglog(Nvals, d1v, 'd-', 'LineWidth', 3,'MarkerSize', 10,'Color', colorp);
   loglog(Nvals, yforplot2, 'LineWidth', 2,'Color', colorg)
text(100,1/10000,'1/N^2', 'Color', colorg, 'FontSize', 16);
hold off
legend('L^{\infty}', 'L^2', 'L^1');
axis([min(Nvals) Nvals(end-1) 10^(-6) 10]);
set(gca, 'FontSize', 12);
xlabel('N', 'FontSize', 18)
ylabel('Error', 'FontSize', 18)
title('v-Error, $\beta = 0.00018$, $\epsilon_1 = Cds^{1/2}$, $\epsilon_2 = C(ds)^{2/5}$, using blob $\psi$',...
    'FontSize', 18,'Interpreter','latex')
