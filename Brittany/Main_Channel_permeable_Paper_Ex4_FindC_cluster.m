
%Example 3 of Cortez, Fluids 2021
%Channel with inflow and part of membrane permeable

%Developed by Ricardo Cortez, Brittany Leathers, and Michaela Kubacki
%July 2024

%Find best C1 for eps=Cds^(1/p) for different values of p
%And for different Flow Regimes

clear all 
% close all

%% Parameters to set 

%setting the viscosity
mu = 1; 
Nvals=10*2.^(0:5);
% Nvals=[10 20 40 80 160 320];
%number of points on boundary where velocity is set and force is computed 

%constant ep=c ds^(1/p)
cvals=1:0.01:4.5;
p=1;

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

usolutions=cell(length(Nvals), length(cvals)); %for storing u solutions to compare
vsolutions=cell(length(Nvals), length(cvals));
%%
for j=1:length(Nvals)
    N=Nvals(j)


for i=1:length(cvals)
   
%discretization of channel
ds = (ymax-ymin)/N;

%%
%regularization parameter
ep=cvals(i)*ds^(1/p);

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


%%
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

% %calculate on grid
% dx=(ymax-ymin)/160;
% % dx=ds;
% % dxs(i)=dx;
% Nx=round(xmax/dx);
% Ny=round(ymax/dx);
% xg=dx*(0:Nx-1)+xmin;
% yg=dx*(0:Ny-1)+ymin;
% [xg,yg]=ndgrid(xg,yg);
% xgv=reshape(xg, Nx*Ny,1);
% ygv=reshape(yg, Nx*Ny,1);
% Ugrid=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[xgv,ygv],ep,mu, blob, ds);
% ug=reshape(Ugrid(:,1), Nx, Ny); 
% vg=reshape(Ugrid(:,2), Nx, Ny); 
% speed=sqrt(ug.^2+vg.^2);
% usolutions{i}=ug;
% vsolutions{i}=vg;
% 
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

error(j,i)=Rin+Rtop+Rout;

end
end
%%
save('error_p1_no.mat', 'error', 'Nvals', 'cvals')