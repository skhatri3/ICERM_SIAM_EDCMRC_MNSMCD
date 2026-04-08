%% Segment Method on the Infinite Channel Problem
% Flow through permeable channel using exact solutions on left and right
% boundaries from Bernardi et al 2023.
%
% Using the segment method with Stokeslets and Source Doublets
%
% Function files needed
%   segment_int_terms.m
%   seg_reg_funcs.m
%   st_segments_velocitytoforce.m
%   st_segments_forcetovelocity.m
%   sd_segments_velocitytogforce.m
%   sd_segments_gforcetovelocity.m
%
% Developed by Michaela Kubacki March 2026, based off of the Infinite
% Channel MRS codes.

%% Parameters to set

% Model Parameters
Da = 0.4; % Darcy Number
mu = 1; % Viscosity
blob_num = 1; % blob choice
C = 0.039; % Coefficient for beta = C*Da*ep, use 0.13 for phi, and 0.039 for psi

% Segment Discretization
Ny1 = 160; % Number of y-points on top/bottom boundary (including corners)
Ny2 = 51;  % Number of y-points on inlet/outlet (don't include corners)

Nseg = 2*(Ny1+Ny2); % Total number of segments on the boundary 

% Evaluation Grid
Nx1 = 200; % Number of target points in x dir for full channel calculation
Nx2 = 80; % Number of target poitns in y dir for full channel calculation

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


% Discretization step
ds_x = (xmax - xmin)/(Ny1-1);
ds_y = (ymax - ymin)/(Ny2+1);

% Regularization parameter -- Pick ep << h for segments
ep =ds_y/10;
beta_value = Da*ep*C;

% Discretize Channel Walls 
% top and bottom
stb = xmin:ds_x:xmax;
stb = stb';
% left and right
slr = ymin+ds_y:ds_y:ymax-ds_y;
slr = slr';

% Define coordinates on each wall (x-coord,y-coord)
y1_top = stb;   y2_top = ymax*ones(size(stb)); % top wall coordinates (y1_top,y2_top)
y1_bot = flip(stb);   y2_bot = ymin*ones(size(stb)); % bottom wall coordinates (y1_bot,y2_bot)
y1_left = xmin*ones(size(slr));   y2_left = slr; % left wall coordinates (y1_left,y2_left)
y1_right = xmax*ones(size(slr));   y2_right = flip(slr); % right wall coordinates (y1_right,y2_right)
y1 = [y1_top; y1_right; y1_bot; y1_left]; %x-coordinates of all boundary points
y2 = [y2_top; y2_right; y2_bot; y2_left]; %y-coordinates of all boundary points

% Indices of permeable region (where the boundary velocity is unknown)
%idx = find(y2 == -1 | y2 == 1); % corners part of permeable boundary
idx = find(-perm_factor*Lx<y1 & y1<perm_factor*Lx); % < since corners are part of solid boundary
%solidI = find(y1 == xmin | y1 == xmax);
beta = zeros(length(y1),1);
beta(idx) = beta_value;


% Define normal vectors
normals_top = zeros(Ny1-1,2); % unit normals for top:
normals_top(:,2) = 1;
normals_bot = zeros(Ny1-1,2); % unit normals for bottom
normals_bot(:,2) = -1;
normals_left = zeros(Ny2+1,2); % unit normals for left  wall 
normals_left(:,1) = -1;
normals_right = zeros(Ny2+1,2); % unit normals for right wall 
normals_right(:,1) = 1;
% normals on full boundary
normals = [normals_top; normals_right; normals_bot; normals_left];

% Points where velocity will be computed within the channel
xx1 = linspace(xmin,xmax,Nx1);
xx2 = linspace(ymin,ymax,Nx2);
dx_g = (xmax-xmin)/Nx1;
dy_g = (ymax - ymin)/Nx2;
[x1m, x2m] = ndgrid(xx1, xx2);
x1 = x1m(:); %x-coords of points where computing velocity
x2 = x2m(:); %y-coords of points where computing velocity
[x1gg, x2gg] = meshgrid(xx1, xx2);
x = [x1,x2];

% Create segment data matrix [left-idx, right-idx, ep, beta-left, beta-right]
% Order of segments: top (l->r), right (t->b), bottom (r->l), left (b->t)
% Must use a continuous ordering of segment points
segData = zeros(Nseg,4);
segData(:,1) = (1:Nseg)';
segData(:,2) = [(2:Nseg)'; 1];
segData(:,3) = ep*ones(Nseg,1); % epsilon value for segment k (using 1 value)
segData(:,4) = beta(segData(:,1)); % beta value at starting points
segData(:,5) = beta(segData(:,2)); % beta value at ending points

% % Segment points (starting points of each segment)
y1seg = y1(segData(:,1)); % x-coords
y2seg = y2(segData(:,1)); % y-coords

beta_start = segData(:,4);
beta_end = segData(:,5);
beta_seg = (beta_start + beta_end)/2;

% % Plotting segments
% figure
% hold on
% axis([-4,4,-1.5,1.5])
% for k = 1:Nseg
%     y = [y1(segData(k,1)) y2(segData(k,1));
%     y1(segData(k,2)) y2(segData(k,2))];
%     plot(y(:,1),y(:,2),'-','linewidth',3)
%     quiver(y1(segData(k,1)),y2(segData(k,1)),normals(k,1),normals(k,2))
%     pause(0.2)
% end
% 
% title('Segments')

%% Stokeslets Only
% Using exact velocities to recover forces using segment method, then
% calculate velocities everywhere using segment method (checking that
% segment method works).

% Using exact velocity from true solution at the segment points
[u1,u2,~] = permeablechannelexact(y1,y2,Da);
% Velocity in segment order
u1seg = u1(segData(:,1));
u2seg = u2(segData(:,1));
useg = [u1seg, u2seg];

% Calculate forces using segment method
[fseg_temp] = st_segments_velocitytoforce([y1,y2],useg,mu,segData,blob_num);
f1seg = fseg_temp(:,1);
f2seg = fseg_temp(:,2);

% Calculate Channel Velocities (checking method)
[U_st_g] = st_segments_forcetovelocity(x,[y1seg,y2seg],fseg_temp,mu,segData,blob_num);

u1_seg_m = reshape(U_st_g(:,1),size(xx1,2),size(xx2,2)); %x-coords of computed velocities
u2_seg_m = reshape(U_st_g(:,2),size(xx1,2),size(xx2,2)); %y-coords of computed velocities
[uexact,vexact,~] = permeablechannelexact(x1gg,x2gg,Da);
u_error_seg = uexact - u1_seg_m';
v_error_seg = vexact - u2_seg_m';

% Make streamline plot
plot_streamline(u1_seg_m,u2_seg_m,x1gg,x2gg,uexact,vexact,y1,y2, 'Stokeslet Seg')

% Errors
[Stokeslet_Seg_errors] = calculate_errors(u_error_seg,v_error_seg,xx1,xx2,ds_x,ds_y,'Stokeslet Seg')


%% Permeable Problem
% Now we run through the algorithm for the permeable channel problem using
% segments.

% STEP 1: Find g using u = utrue on solid bd and u = 0 on permeable bd
% Using exact velocity from true solution at the segment points
[u1,u2,~] = permeablechannelexact(y1,y2,Da);
u = [u1, u2];
% Velocity in segment order
u1seg = u1(segData(:,1));
u2seg = u2(segData(:,1));

% Make velocity with 0s in permeable sections
u1b = u1; %x-coordinates of all boundary velocities
u2b = u2; %y-coordinates of all boundary velocities
u1b(idx) = 0; % set boundary velocities on permeable section to zero
u2b(idx) = 0;

% THIS IS WHERE THE PROBLEM IS!!!!
[gseg] = sd_segments_velocitytogforce([y1,y2],[u1b,u2b],normals,mu,segData,blob_num);
%[gseg] = st_segments_velocitytoforce([y1,y2],[u1b,u2b],mu,segData,blob_num); % why does this give better results?
g1seg = gseg(:,1);
g2seg = gseg(:,2);

% STEP 2: Recover boundary velocity on permeable boundary using g
y1b = y1(idx); % boundary points of permeable section of boundary
y2b = y2(idx);

% Recover missing boundary velocity pieces on permeable section
[useg_perm] = sd_segments_gforcetovelocity([y1b,y2b],[y1,y2],gseg,normals,mu,segData,blob_num);
%[useg_perm] = st_segments_forcetovelocity([y1b,y2b],[y1,y2],gseg,mu,segData);
useg_perm1 = useg_perm(:,1);
useg_perm2 = useg_perm(:,2);

u1b_seg = u1b;
u2b_seg = u2b;

% Replace the 0s with the newly computed velocities in permeable region
u1b_seg(idx) = useg_perm1;
u2b_seg(idx) = useg_perm2;

% STEP 3: Use boundary velocities to solve for forces on boundary
ub_seg = [u1b_seg, u2b_seg];
[fseg] = st_segments_velocitytoforce([y1,y2],ub_seg,mu,segData,blob_num);

% STEP 4: Use forces on boundary to calculate velocities throughout channel
[useg_g] = st_segments_forcetovelocity(x,[y1seg,y2seg],fseg,mu,segData,blob_num);

%% Figures


% Source doublets (g)
figure;
scalefactor = 0.001;
sk = 2;
quiver(y1(1:sk:end),y2(1:sk:end),g1seg(1:sk:end)*scalefactor,g2seg(1:sk:end)*scalefactor,'AutoScale','off')
axis equal;
axis([-Lx-1,Lx+1,-1.5,1.5])
title('Distribution of Source Doublets (g)')

% Exact vs recovered boundary velocity
sk = 1;
scalefactor = 0.05; % for the boundary velocity figures
figure
subplot(2,1,1);
hold on;
plot(y1,y2,'k.')
hold on
quiver(y1seg(1:sk:end),y2seg(1:sk:end),u1b_seg(1:sk:end)*scalefactor,u2b_seg(1:sk:end)*scalefactor,'b','AutoScale','off')
title('Recovered Bd Velocities on Boundaries using segment method')
set(gca,'Fontsize',14)
hold off;
axis equal;
axis([-Lx-1,Lx+1,-1.5,1.5])
subplot(2,1,2);
plot(y1,y2,'k.')
hold on
quiver(y1(1:sk:end),y2(1:sk:end),u1(1:sk:end)*scalefactor, u2(1:sk:end)*scalefactor,'r','Autoscale','off')
hold off
title('Exact Boundary Velocities')
set(gca,'Fontsize',14)
axis equal
axis([-Lx-1,Lx+1,-1.5,1.5])

% Permeable Segment Errors: Table and Plots
u1_seg_m = reshape(useg_g(:,1),size(xx1,2),size(xx2,2)); %x-coords of computed velocities
u2_seg_m = reshape(useg_g(:,2),size(xx1,2),size(xx2,2)); %y-coords of computed velocities
[uexact,vexact,~] = permeablechannelexact(x1gg,x2gg,Da);
u_error_seg = uexact - u1_seg_m';
v_error_seg = vexact - u2_seg_m';

[Permeable_Seg_errors] = calculate_errors(u_error_seg,v_error_seg,xx1,xx2,ds_x,ds_y,'Permeable Seg')

plot_streamline(u1_seg_m,u2_seg_m,x1gg,x2gg,uexact,vexact,y1,y2,'Permeable Seg')

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
%% Errors
function [errors,Emax,EL1,EL2] = calculate_errors(u_error,v_error,xx1,xx2,dx,dy, method_name)

u_error_max = max(max(abs(u_error)));
v_error_max = max(max(abs(v_error)));
u_error_L1 = dx*dy*sum(sum(abs(u_error)));
v_error_L1 = dx*dy*sum(sum(abs(v_error)));
u_error_L2 = sqrt(dx*dy*sum(sum(u_error.^2)));
v_error_L2 = sqrt(dx*dy*sum(sum(v_error.^2)));

ColNames = {'Max','L1','L2'};
Emax = [u_error_max; v_error_max];
EL1 = [u_error_L1; v_error_L1];
EL2 = [u_error_L2; v_error_L2];
errors = table(Emax, EL1, EL2,'VariableNames',ColNames);

figure
subplot(2,1,1);
surf(xx1, xx2, u_error, 'EdgeColor', 'none');
colorbar;
title(['Error in u (horizontal) velocity -- ' method_name]);
xlabel('x');
ylabel('y');
axis equal;
view(2); % View from above for 2D representation

subplot(2,1,2);
surf(xx1, xx2, v_error, 'EdgeColor', 'none');
colorbar;
title(['Error in v (vertical) velocity -- ' method_name]);
xlabel('x');
ylabel('y');
axis equal;
view(2); % View from above for 2D representation

end
%% Streamline plots
function plot_streamline(u1_seg_m,u2_seg_m,x1gg,x2gg,uexact,vexact,y1,y2,method_name)

figure
plot(y1, y2, 'k.')
hold on

% Computed streamlines
hh1_comp = streamline(x1gg, x2gg, u1_seg_m', u2_seg_m', x1gg(1:end, 1), x2gg(1:end, 1)); % streamlines starting at left wall
hh2_comp = streamline(x1gg, x2gg, -u1_seg_m', -u2_seg_m', x1gg(1:end, end), x2gg(1:end, end)); % streamlines starting at right wall
% Exact streamlines
hh1_exact = streamline(x1gg, x2gg, uexact, vexact, x1gg(1:end, 1), x2gg(1:end, 1)); % streamlines starting at left wall
hh2_exact = streamline(x1gg, x2gg, -uexact, -vexact, x1gg(1:end, end), x2gg(1:end, end)); % streamlines starting at right wall

strmLW = 2;
RGB = orderedcolors("gem");
set(hh1_comp, 'Color', RGB(1,:), 'linewidth', strmLW);
set(hh2_comp, 'Color', RGB(1,:), 'linewidth', strmLW);
set(hh1_exact, 'Color', RGB(3,:), 'LineStyle','--', 'linewidth', strmLW);
set(hh2_exact, 'Color', RGB(3,:), 'LineStyle','--', 'linewidth', strmLW);

title(['Computed (blue) and exact (yellow) streamlines - ' method_name])
set(gca, 'FontSize', 14)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
axis equal
axis([-4,4,-1.5,1.5])
end
