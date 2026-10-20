%Example 3 of Cortez, Fluids 2021
%Channel with inflow and part of membrane permeable

%Developed by Ricardo Cortez, Brittany Leathers, and Michaela Kubacki
%July 2024

% This script is for determining the relationship between beta and epsilon
% so that we can adjust beta appropriately to preserve the flow when
% changing epsilon.

clear all
close all

addpath('./stokeslets_codes');

%% Parameters to set

%setting the viscosity
mu = 1;

%number of points on boundary where velocity is set and force is computed
N = 160;

%blob choice
blob=2;

%permeability coefficient (assuming constant for the permeable region)
%b=0.0001;

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
%regularization parameter
ep1 = 0.010;
ep2 =ep1;%0.05/sqrt(5);%1*ds;


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
% y's for permeable region
y1b=y1(I);
y2b=y2(I);
% Target Coordinates for outlet
x2_out= (ymin+ds/2:ds:ymax-ds/2)';
x1_out= xmax*ones(size(x2_out));
%unit normals
normals=[normals_top; normals_bot; normals_side];

% Set beta vector in for-loop
%beta vector for function;
%beta= zeros(length(y1),1);
%beta(I)=b;

%% Initial Values for bisection

b1 = 0.00007;
b2 = 0.00009;

Rtarget = 0.009889411831155; % blob 2 for ep = 0.05/sqrt(5)

Routb1 =0.089610154420317;

Routb2 = -0.074588847964431;

% For b = 0.0001
%Rtarget = 0.302653177713578; % for ep = 0.05/sqrt(5), b = 0.0001, blob 2
%Rtarget = -0.210667757702758; % for ep = 0.05/sqrt(5), b = 0.00018, blob 2
%Routb1 = -0.166415298257469;
%Routb2 = -0.308303496141013;

% For b = 0.0018
%Rtarget = 0.009889411831151; % b = 0.00018, ep = 0.05/sqrt(5), blob 2;




Rdiff1 =  Routb1-Rtarget % for current ep and b1
Rdiff2 = Routb2-Rtarget % for current ep and b2
% check these have opposite signs

l = b2-b1;
tol = 1e-6;
itMax = 30;
it_count = 0;
for k = 1:itMax
    it_count = it_count + 1;
    b = (b1+b2)/2;
    l = (b2-b1)/2;
    %%
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

    %beta vector for function using current b values;
    beta= zeros(length(y1),1);
    beta(I)=b;

    %% Compute solution for current beta
    %compute g
    g=RegStokeslets2D_velocityto_gforce_permeable_diff_ep([y1,y2],[y1,y2],...
        [u1,u2],ep1,ep2,mu, blob, I, beta, normals);

    %Find velocity in permeable region:

    [u_beta] = RegStokeslets2D_gtovelocity([y1,y2],g, [y1b,y2b],ep2,mu,blob, beta, normals);

    %Put in the new velocities for the permeable part
    %permeable part : temporary velocity for finding g
    u1(I)=u_beta(:,1);
    u2(I)=u_beta(:,2);

    %computing the force
    f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep1,mu, blob);
    f1 = f(:,1);
    f2 = f(:,2);

    %outlet
    uu = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],...
        [x1_out,x2_out],ep1,mu, blob);
    Rout=-ds*sum(dot(normals_side, uu))

    Rdiff = Rout-Rtarget % current difference in outflows

    current = [it_count, b1, b2, b, Rdiff, l];
    disp(current)

    if (abs(Rdiff)< tol)
        disp('Bisection Converged')
        break
    end
    if sign(Rdiff)~=sign(Rdiff1) % root is in [b1 b]
        b2 = b;
        Rdiff2 = Rdiff;
    else % root is in [b b2]
        b1 = b;
        Rdiff1 = Rdiff;
    end
    if (it_count >= itMax)
        disp('Did not converge within tol. Exceeded itMax')
    end
end
