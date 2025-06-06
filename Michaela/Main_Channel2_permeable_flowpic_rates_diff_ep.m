%Example 3 of Cortez, Fluids 2021
%Channel with inflow and part of membrane permeable

%Developed by Ricardo Cortez, Brittany Leathers, and Michaela Kubacki
%July 2024


clear all 
close all

format long

currDate = datestr(datetime);

mkdir(currDate)
currPath = [pwd,'/',currDate];

resultsFile = fopen([currPath,'/results.txt'],'w');
fprintf(resultsFile,['Results for ',currDate]);
fprintf(resultsFile,'\nParameters\n');


%% Parameters to set 

%setting the viscosity
mu = 1; 

%number of points on boundary where velocity is set and force is computed  
N = 160;
% Nvalues = 10*2.^(0:4);

%blob choice
blob=2;

% C-values for blob 2
C1 = 3.5; % BL C1-value
C2 = 0.06; % BL C2-value
%C1 = 3.45; % MK C1-value
%C2 = 0.0123; % MK C2-value

fprintf(resultsFile,['mu = ',num2str(mu),', using blob ',num2str(blob)]);
fprintf(resultsFile,['\nN = ',num2str(N)]);

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

%discretization of channel
ds = (ymax-ymin)/N;
%%
%regularization parameters
%ep1 =C1*ds;%1*ds;
ep1 = C1*ds;
ep2 = C2*(ds)^(1/3);

%permeability coefficient (assuming constant for the permeable region) 
% for blob 1 (phi)
%b = 0.011843080738230*ep2; % for b = 0.0001 (beta 1)
%b = 0.004247164369662*ep2; % for b = 0.00018
%b = 0.0282306994128518*ep2; % for b = 0.00024 flow


%permeability coefficient (assuming constant for the permeable region) 
% for blob 2 (psi)
% b=0.010733*ep2; % for b = 0.00024 (beta 3)
 b = 0.00805*ep2; % for b = 0.00018 (beta 2)
% b = 0.004472*ep2; % for b = 0.0001 (beta 1)
fprintf(resultsFile,[', \beta =  ', num2str(b)]);
fprintf(resultsFile,[', epsilon1 = C1ds = ', num2str(ep1), ', epsilon2 = C2(ds)^(1/3)=', num2str(ep2)]);

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
g=RegStokeslets2D_velocityto_gforce_permeable_diff_ep([y1,y2],[y1,y2],...
    [u1,u2],ep1,ep2,mu, blob, I, beta, normals);

%Find velocity in permeable region:
y1b=y1(I);
y2b=y2(I);
[u_beta] = RegStokeslets2D_gtovelocity([y1,y2],g, [y1b,y2b],ep2,mu,blob, beta, normals);
% 
% u1=u1+u_beta(:,1);
% u2=u2+u_beta(:,2);

%Put in the new velocities for the permeable part
%permeable part : temporary velocity for finding g
u1(I)=u_beta(:,1);
u2(I)=u_beta(:,2);

%computing the force 
f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep1,mu, blob);
f1 = f(:,1);
f2 = f(:,2);  

% Calculate velocities on right hand side too
Ufull=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[y1f,y2f],ep1,mu, blob, ds);
ufull=Ufull(:,1); vfull=Ufull(:,2);

%calculate on grid
[xgg,ygg] = meshgrid(ds:4*ds:Lx-ds, ds:4*ds:Ly-ds); xg=xgg(:); yg=ygg(:);
Ugrid=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[xg,yg],ep1,mu, blob, ds);
ug=reshape(Ugrid(:,1), length(ds:4*ds:Ly-ds), length(ds:4*ds:Lx-ds)); 
vg=reshape(Ugrid(:,2), length(ds:4*ds:Ly-ds), length(ds:4*ds:Lx-ds)); 
speed=sqrt(ug.^2+vg.^2);
%% Plot figure
sk=2;
figure;%subplot(211)
plot(y1f,y2f,'k.'),hold on
% quiver(xgg,ygg,ug,vg,0.8,'r','LineWidth',1)
quiver(y1f(1:6*sk:end),y2f(1:6*sk:end),ufull(1:6*sk:end),vfull(1:6*sk:end),0,'k','LineWidth',2)
% quiver(xe(1:8:end),ye(1:8:end),usuck(1:8:end),vsuck(1:8:end),0,'r')
surf(xgg,ygg,-speed),view(2),shading interp
hh1=streamline(xgg,ygg,ug ,vg ,xgg(1:end,1 ),ygg( 1:end,1 ));
hh2=streamline(xgg,ygg,ug ,vg ,xgg(1:end,end ),ygg( 1:end,end ));
set(hh1,'Color','black');hold off,axis equal,
set(hh2,'Color','black');hold off,axis equal,axis([-0.10 5.5 -0.55 1.75  ])
title(['\epsilon_1 = ',num2str(ep1), ', \epsilon_2 = ',num2str(ep2), ' (\beta = ',num2str(b), ')'])
% clim([0 1]);
colorbar('Ticks',[-1 , -0.8, -0.6,-0.4,-0.2,0 ],...
    'TickLabels',{'1','0.8','0.6','0.4','0.2','0'},...
    'Direction','reverse')
hold off  
set(gca, 'FontSize', 16)
fig_title = ['epsilon_1 = ',num2str(ep1),', \epsilon_2 = ',num2str(ep2),'\beta = ',num2str(b),'.fig'];
savefig([currPath,'/',fig_title])

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
    [x1_out,x2_out],ep1,mu, blob);
Rout=-ds*sum(dot(normals_side, uu));

NetFlow = Rin+Rtop+Rout;

fprintf(resultsFile,'\n Flow Rates\n');
fprintf(resultsFile,['Rinlet = ',num2str(Rin),', Rtop = ',num2str(Rtop),', Routlet = ',num2str(Rout)]);
fprintf(resultsFile,['\n NetFlow = ',num2str(NetFlow)]);
