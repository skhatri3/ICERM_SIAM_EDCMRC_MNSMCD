% Flow through permeable channel using exact solutions on left and right
% boundaries from Bernardi et al 2023.

% Validation of segment method on 1 segment of the boundary, comparing
% results to what we obtain on the segment using MRS with blob phi. 
%
% Developed by Michaela Kubacki and Ricardo Cortez March 2026


clear all
close all

%% Parameters to set

% Model Parameters
Da = 0.4; % Darcy Number
mu = 1; % Viscosity
blob_num = 2; % blob choice

% Number of segment and target points
%N = 160; % Number of source points (along top  and bottom boundaries)

Ny1 = 30; % Number of y-points on top/bottom boundary
% Ny1 + 1 segments on top/bottom
Ny2 = 10;  % Number of y-points on inlet/outlet (includes corners)
% Ny2 - 1 segments on left/right
Nseg = 2*(Ny1+Ny2); % Total number of segments on the boundary

Nx1 = 40*5; % Number of target points in x direction for full channel calc
Nx2 = 80; % Number of target points in y direction for full channel calc

% Setting up channel

%            1 - top (perm)
%         +""""""""""""""""""""""""""""""+
%         |                              |
%  3-left |                              | 4-right
% (known) |                              | (known)
%         |                              |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |   
%         +""""""""""""""""""""""""""""""+
%         2 - bottom (perm)

% Channel Geometry
Lx = pi;
Ly = 1;
xmin = -Lx;
xmax = Lx;
ymin = -Ly;
ymax = Ly;
perm_factor=1;  % will make permeable section  -factor*Lx < x < factor*Lx
scalefactor = 0.05; % for the boundary velocity figures

% Discretization step
ds_x = (xmax - xmin)/(Ny1+1);
ds_y = (ymax - ymin)/(Ny2-1);


% Regularization parameter
ep = ds_y/10;%ds_y;
beta_value = Da*ep*0.04;
%beta_value = @(x) -0.6 - 0.35*(((pi-x).*(pi+x))/pi^2).^(4/6);
% Discretize Channel Walls (corners belong to left and right)
% top and bottom
stb = xmin+ds_x:ds_x:xmax-ds_x;
stb = stb';
% left and right
slr = ymin:ds_y:ymax;
slr = slr';

% Define coordinates on each wall (x-coord,y-coord)
y1_top = stb;   y2_top = ymax*ones(size(stb)); % top wall coordinates (y1_top,y2_top)
y1_bot = stb;   y2_bot = ymin*ones(size(stb)); % bottom wall coordinates (y1_bot,y2_bot)
y1_left = xmin*ones(size(slr));   y2_left = slr; % left wall coordinates (y1_left,y2_left)
y1_right = xmax*ones(size(slr));   y2_right = slr; % right wall coordinates (y1_right,y2_right)
y1 = [y1_top; y1_bot; y1_left; y1_right]; %x-coordinates of all boundary points
y2 = [y2_top; y2_bot; y2_left; y2_right]; %y-coordinates of all boundary points

% Indices of permeable region (where the boundary velocity is unknown)
%idx = find(-perm_factor*Lx<y1 & y1<perm_factor*Lx);
%solidI = find(y1 == xmin | y1 == xmax);
%beta = zeros(length(y1),1);
beta = zeros(Nseg,1);
idx = 1:(2*Ny1+2);
beta(idx) = beta_value;

segData = zeros(Nseg,4);

% x-coord
segData(:,1) = [2*Ny1+Ny2; (1:Ny1)'; 2*Ny1+1; (Ny1+1:2*Ny1)'; 
    (2*Ny1+1:2*Ny1+Ny2-1)'; (2*Ny1+Ny2+1:2*Ny1+2*Ny2-1)'];
% y-coord 
segData(:,2) = [(1:Ny1)'; 2*Ny1+2*Ny2; (Ny1+1:2*Ny1)'; 2*Ny1+Ny2+1; 
    (2*Ny1+2:2*Ny1+Ny2)'; (2*Ny1+Ny2+2:2*Ny1+2*Ny2)'];
segData(:,3) = ep*ones(Nseg,1);
segData(:,4) = beta.*ones(Nseg,1);

% Points where velocity will be computed within the channel
xx1 = linspace(xmin+ds_x/2,xmax-ds_x/2,Nx1);
xx2 = linspace(ymin+ds_y/2,ymax-ds_y/2,Nx2);
[x1m, x2m] = ndgrid(xx1, xx2);
x1 = x1m(:); %x-coords of points where computing velocity
x2 = x2m(:); %y-coords of points where computing velocity

% Stokeslet and Source Doublet segment velocities using just one segment,
% with multiple evaluation points.

I = eye(2,2);
k = 1; % Specify segment
ep = segData(k,3); % Epsilon value on segment
% y = segment starting and ending points
y = [y1(segData(k:k+1,1)) y2(segData(k:k+1,1))];
% Define evaluation points
x1 = linspace(-pi,pi,10)';
x2 = 0.5*ones(size(x1));
x = [x1 x2];


%[Q0,Q1,T02,T12,T22,T32] = st_segment_terms(x,y,ep);
% Specify forces at endpoints of segment
fj = [0;1];
fk = [0;2];
u_st_seg = zeros(size(x));
u_sd_seg = zeros(size(x));
for m = 1:length(x)
    % Calculate Stokeslet segment velocity from forces
    [Q0,Q1,T02,T12,T22,T32] = st_segment_terms(x(m,:),y,ep);
    x0 = (x(m,:) - y(1,:))';
    v = ((y(1,:)-y(2,:)))';
    L = norm(v);
    M2 = ((ep^2*T12 - Q1)*I + T12.*(x0*x0') + T22.*(x0*v' + v*x0') + T32.*(v*v'));
    M1 = ((ep^2*T02 - Q0)*I + T02.*(x0*x0')+ T12.*(x0*v' + v*x0') + T22.*(v*v')) - M2;
    % Calculate Source Doublet segment velocity from forces
    u_st_seg(m,:) = ((M1*fj + M2*fk)/(4*pi))';
    [T02,T04,T06,T12,T14,T16,T24,T26] = sd_segment_terms(x(m,:),y,ep);
    n = -[-v(2,1); v(1,1)]/L; % Outward Unit normal vector for this segment
    nx0 = n'*x0; 
    x0n = x0*n';
    vn = v*n';
    nn = n*n';
    D2 = ((T12 + ep^2*T14).*nn + nx0.*(-2*T14 -4*ep^2*T16).*x0n + nx0.*(- 2*T24 - 4*ep^2*T26).*vn);
    D1 = (T02 + ep^2*T04).*nn + nx0.*(-2*T04 - 4*ep^2*T06).*x0n + nx0.*(-2*T14 - 4*ep^2*T16).*vn  - D2;
    u_sd_seg(m,:) = -((D1*fj + D2*fk)/(2*pi));
end

% Check with MRS using a linear distribution of force densities
h = 1/32;
yc = (y(1,:)' - (h/2:h:1-h/2).*v)';
f = (fj + (h/2:h:1-h/2).*(fk-fj))';
wt = [h*ones(size(yc))];
beta = ones(size(yc))';
normals = zeros(length(yc),2);
% Segment is on top boundary, normals are [0,1]
normals(:,1) = zeros(length(yc),1);
normals(:,2) = ones(length(yc),1);
% Calculate velocities from given forces
[u_mrs_check] = RegStokeslets2D_forcetovelocity(yc,f,x,ep,mu,1,wt);
[u_sd_check] = RegStokeslets2D_gtovelocity(yc,f,x,ep,mu,1, beta, normals, wt);

StComp = [u_st_seg-u_mrs_check]
SDcomp = [u_sd_seg - u_sd_check]

%% Exact solutions
function [ug,vg,pg] = permeablechannelexact(x,y,Da)

% Newton's method to find lambda1
L1=sqrt(2*Da);
f=@(L) tan(L)^2-tan(L)/L+1-2*Da;
L1 = fzero(f,L1);
%disp(['lambda_1 = ',num2str(L1)])

% Compute g and g' from FB paper
g  = cot(L1)*(Da-0.5)   *sin(L1*y)+0.5*y.*cos(L1*y);
gp = cot(L1)*(Da-0.5)*L1*cos(L1*y)+0.5*cos(L1*y)-0.5*L1*y.*sin(L1*y);

% Compute exact solution using g and g'
C = -4;
pg = C*sinh(L1*x) .* cos(L1*y);
ug = -C/L1 * cosh(L1*x) .* gp;
vg =  C    * sinh(L1*x) .* g ;

end