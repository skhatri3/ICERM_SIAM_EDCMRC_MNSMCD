%Example 4 of Cortez, Fluids 2021
%Channel with inflow and part of membrane permeable

%Developed by Ricardo Cortez, Brittany Leathers, and Michaela Kubacki
%July 2024

%Refinement studies as refine discretization of channel,
%letting Beta=A eps and eps=C ds^(1/p)
%Refinement study on net flow out 
%Refinement study on L1,2,inf norms of sequential differences in velocity
%   solutions calculated on a grid (grid constant, not refined) in the channel


clear all 
% close all

%% Parameters to set 

%setting the viscosity
mu = 1; 

%number of points on boundary where velocity is set and force is computed 
Nvals = 20*2.^(-1:4);    

%blob choice
blob=2;

%constant determining inflow velocity profile
a=4;

  

%Degree to for convergence (eps=C ds^(1/p))
p=4;
 
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


%%
for i=1:length(Nvals)
    N=Nvals(i)

%discretization of channel
ds = (ymax-ymin)/N;
ds_s(i)=ds;
%%
%regularization parameter
C=0.05/sqrt(5)*160^(1/p);
ep =C*ds^(1/p); 

%permeability coefficient (assuming constant for the permeable region)
b=0.00018*sqrt(5)/0.05*ep;

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

%coordiantes of boundary points
y1 = [y1_top;y1_bot; y1_side];
y2 = [y2_top;y2_bot; y2_side];
%For use at very end:
y2right = (ymin+ds/2:ds:ymax-ds/2)';
y1right = xmax*ones(size(y2_side));
y1f=[y1; y1right];
y2f=[y2; y2right];
%indices of permeable region
I=find(y1_top<2/3*Lx & y1_top>1/3*Lx);
I2=find(y1_top>=2/3*Lx | y1_top<=1/3*Lx);
%unit normals
normals=[normals_top; normals_bot; normals_side];
%beta vector for function;
beta= zeros(length(y1),1); 
beta(I)=b;



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
%compute g
g=RegStokeslets2D_velocityto_gforce_permeable([y1,y2],[y1,y2],...
    [u1,u2],ep,mu, blob, I, beta, normals);

%Find velocity in permeable region:
y1b=y1(I);
y2b=y2(I);
[u_beta]=RegStokeslets2D_gtovelocity([y1,y2],g, [y1b,y2b],...
    ep,mu, blob, beta, normals);
% 
% u1=u1+u_beta(:,1);
% u2=u2+u_beta(:,2);

%Put in the new velocities for the permeable part
%permeable part : temporary velocity for finding g
u1(I)=u_beta(:,1);
u2(I)=u_beta(:,2);

%computing the force 
f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep,mu, blob);
f1 = f(:,1);
f2 = f(:,2);  

% Calculate velocities on right hand side too
Ufull=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[y1f,y2f],ep,mu, blob, ds);
ufull=Ufull(:,1); vfull=Ufull(:,2);

%calculate on grid
dx=(ymax-ymin)/160;
% dx=ds;
% dxs(i)=dx;
Nx=round(xmax/dx);
Ny=round(ymax/dx);
xg=dx*(0:Nx-1)+xmin;
yg=dx*(0:Ny-1)+ymin;
[xg,yg]=ndgrid(xg,yg);
xgv=reshape(xg, Nx*Ny,1);
ygv=reshape(yg, Nx*Ny,1);
Ugrid=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[xgv,ygv],ep,mu, blob, ds);
ug=reshape(Ugrid(:,1), Nx, Ny); 
vg=reshape(Ugrid(:,2), Nx, Ny); 
speed=sqrt(ug.^2+vg.^2);
usolutions{i}=ug;
vsolutions{i}=vg;

% xgo=xg';
% ygo=yg';
% % %% Plot figure
% sk=2;
% figure;%subplot(211)
% plot(y1f,y2f,'k.'),hold on
% % quiver(xgg,ygg,ug,vg,0.8,'r','LineWidth',1)
% quiver(y1f(1:6*sk:end),y2f(1:6*sk:end),ufull(1:6*sk:end),vfull(1:6*sk:end),0,'k','LineWidth',2)
% % quiver(xe(1:8:end),ye(1:8:end),usuck(1:8:end),vsuck(1:8:end),0,'r')
% surf(xg,yg,-speed),view(2),shading interp
% hh1=streamline(xg',yg',ug' ,vg' ,xgo(1:3*sk:end,1 ),ygo( 1:3*sk:end,1 ));
% hh2=streamline(xg',yg',ug' ,vg' ,xgo(1:3*sk:end,end ),ygo( 1:3*sk:end,end ));
% set(hh1,'Color','black');hold off,axis equal,
% set(hh2,'Color','black');hold off,axis equal,axis([-0.10 5.5 -0.55 1.75  ])
% title(['\beta = ',num2str(b)])
% % clim([0 1]);
% colorbar('Ticks',[-1 , -0.8, -0.6,-0.4,-0.2,0 ],...
%     'TickLabels',{'1','0.8','0.6','0.4','0.2','0'},...
%     'Direction','reverse')
% hold off  
% set(gca, 'FontSize', 16)


% Check flow rates 
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
    [x1_out,x2_out],ep,mu, blob);
Rout=-ds*sum(dot(normals_side, uu));

error(i)=Rin+Rtop+Rout;

end
%%

colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

forploty0=0.1;
forplotdx0=20;
dxforplot=1./Nvals;
yforplot=forploty0/forplotdx0^1*dxforplot.^1;
forploty02=20;
forplotdx02=10;
dxforplot2=1./Nvals;
yforplot2=forploty02/forplotdx02^2*dxforplot2.^2;
   
%Refinement study plots
figure;
loglog(Nvals(2:5), abs(error(2:5)),'o-', 'LineWidth', 2.5, 'MarkerSize', ...
       10, 'Color', colorlb)
hold on
% loglog(Nyvals, e2u, 's-', 'LineWidth', 2.5,'MarkerSize', 10, 'Color', colordb);
% loglog(Nyvals, e1u, 'd-', 'LineWidth', 2.5,'MarkerSize', 10,'Color', colorp);
   loglog(Nvals(4:5), yforplot(4:5), 'LineWidth', 2,'Color', colorg)
text(250,1.5/50000,'$1/N$', 'Color', colorg, 'FontSize', 14, ...
    'Interpreter', 'latex');
loglog(Nvals(4:5), yforplot2(4:5), 'LineWidth', 2,'Color', colorg)
text(200,1.6/1000000,'$1/N^2$', 'Color', colorg, 'FontSize', 14, ...
    'Interpreter', 'latex');
hold off
%legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
axis([20 550 5*10^(-7) 10^(-1) ]);
set(gca, 'FontSize', 17);
% legend('$L^{\infty}$', '$L^2$', '$L^1$','Interpreter','latex', 'FontSize', 20);
xlabel('$N$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Net Flow Refinement Study, $\epsilon=3.584*ds$',...
    'FontSize', 18,'Interpreter','latex')  
% yticks([10^(-5) 10^(-3) 10^(0)])


%%
d1u=zeros(1, length(Nvals));
d2u=zeros(1, length(Nvals));
dmaxu=zeros(1, length(Nvals));

d1v=zeros(1, length(Nvals));
d2v=zeros(1, length(Nvals));
dmaxv=zeros(1, length(Nvals));

for i=1:length(Nvals)-1
    %from fine to coarse 
    ufine=usolutions{i+1}; 
    ucoarse=usolutions{i};
    diffu=(ucoarse-ufine);
    % ufinerest=ufine(1:2:end, 1:2:end);
    % diffu=(ucoarse-ufinerest);


    vfine=vsolutions{i+1}; 
    vcoarse=vsolutions{i};
    diffv=(vcoarse-vfine);
    % vfinerest=vfine(1:2:end, 1:2:end);
    % diffv=(vcoarse-vfinerest);
    % dx=dxs(i);

    d1u(1, i)= dx^2*sum(sum(abs(diffu)));
    d2u(1, i)= sqrt(dx^2*sum(sum(diffu(:,:).^2)));
    dmaxu(1, i)= max(max(abs(diffu)));

    d1v(1, i)= dx^2*sum(sum(abs(diffv)));
    d2v(1, i)= sqrt(dx^2*sum(sum(diffv(:,:).^2)));
    dmaxv(1, i)= max(max(abs(diffv)));

end

%%

%% Allnorms refinement study
   
colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

forploty0=1.5;
forplotdx0=10;
dxforplot=1./Nvals;
yforplot=forploty0/forplotdx0^1*dxforplot.^1;
forploty02=50;
forplotdx02=10;
dxforplot2=1./Nvals;
yforplot2=forploty02/forplotdx02^2*dxforplot2.^2;
   
%Refinement study plots
figure;
loglog(Nvals, dmaxu,'o-', 'LineWidth', 3, 'MarkerSize', ...
       10, 'Color', colorlb)
hold on
loglog(Nvals, d2u, 's-', 'LineWidth', 3,'MarkerSize', 10, 'Color', colordb);
loglog(Nvals, d1u, 'd-', 'LineWidth', 3,'MarkerSize', 10,'Color', colorp);
   loglog(Nvals(3:5), yforplot(3:5), 'LineWidth', 2,'Color', colorg)
text(360,1/7000,'1/N', 'Color', colorg, 'FontSize', 16);
% loglog(Nvals(3:5), yforplot2(3:5), 'LineWidth', 1.5,'Color', colorg)
% text(100,1.4/100000,'1/N^2', 'Color', colorg, 'FontSize', 12);
hold off
%legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
legend('L^{\infty}', 'L^2', 'L^1');
axis([10 360 10^(-4) 10^(0)]);
set(gca, 'FontSize', 12);
xlabel('N_x', 'FontSize', 18)
ylabel('Difference Norms', 'FontSize', 18)
title('Horizontal Error Refinement Study',...
    'FontSize', 18,'Interpreter','latex')
   


figure;
   loglog(Nvals, dmaxv,'o-', 'LineWidth', 3, 'MarkerSize', ...
       10, 'Color', colorlb)
hold on
loglog(Nvals, d2v, 's-', 'LineWidth', 3,'MarkerSize', 10, 'Color', colordb);
loglog(Nvals, d1u, 'd-', 'LineWidth', 3,'MarkerSize', 10,'Color', colorp);
   loglog(Nvals(3:5), yforplot(3:5), 'LineWidth', 2,'Color', colorg)
text(360,1/7000,'1/N', 'Color', colorg, 'FontSize', 16);
% loglog(Nvals(3:5), yforplot2(3:5), 'LineWidth', 1.5,'Color', colorg)
% text(100,1.4/100000,'1/N^2', 'Color', colorg, 'FontSize', 12);
hold off
%legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
legend('L^{\infty}', 'L^2', 'L^1');
axis([10 360 10^(-4) 10^(0)]);
set(gca, 'FontSize', 12);
xlabel('N_x', 'FontSize', 18)
ylabel('Difference Norms', 'FontSize', 18)
title('Vertical Error Refinement Study',...
    'FontSize', 18,'Interpreter','latex')


