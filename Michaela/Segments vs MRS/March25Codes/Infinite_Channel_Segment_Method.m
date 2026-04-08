%% Segment Method
% Flow through permeable channel using exact solutions on left and right
% boundaries from Bernardi et al 2023.
%
% Using the segment method with Stokeslets and Source Doublets
%
% Function files needed
%   st_segment_terms.m
%   sd_segment_terms.m
%   st_segments_velocitytoforce.m
%   st_segments_forcetovelocity.m
%   sd_segments_velocitytogforce.m
%   sd_segments_gforcetovelocity
%
% Developed by Michaela Kubacki March 2026, based off of the Infinite
% Channel MRS codes.

%% Parameters to set

% Model Parameters
Da = 0.4; % Darcy Number
mu = 1; % Viscosity
blob_num = 2; % blob choice

Ny1 = 5; % Number of y-points on top/bottom boundary (including corners)
Ny2 = 2;  % Number of y-points on inlet/outlet (don't include corners)

Nseg = 2*(Ny1+Ny2); % Total number of segments on the boundary 

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
scalefactor = 0.05; % for the boundary velocity figures

% Discretization step
ds_x = (xmax - xmin)/(Ny1-1);
ds_y = (ymax - ymin)/(Ny2+1);

% Regularization parameter -- Pick ep << h for segments
ep =ds_y/100;
C = 0.15;
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
idx = find(-perm_factor*Lx<y1 & y1<perm_factor*Lx); % corners part of solid boundary
%solidI = find(y1 == xmin | y1 == xmax);
beta = zeros(length(y1),1);
beta(idx) = beta_value;

% Define normal vectors
normals_top = zeros(length(stb),2); % unit normals for top:
normals_top(:,2) = 1;
normals_bot = zeros(length(stb),2); % unit normals for bottom
normals_bot(:,2) = -1;
normals_left = zeros(length(slr),2); % unit normals for left  wall 
normals_left(:,1) = -1;
normals_right = zeros(length(slr),2); % unit normals for right wall 
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

% % Plotting segments
% figure
% hold on
% axis([-4,4,-1.5,1.5])
% for k = 1:Nseg
%     y = [y1(segData(k,1)) y2(segData(k,1));
%     y1(segData(k,2)) y2(segData(k,2))];
%     plot(y(:,1),y(:,2),'-','linewidth',3)
%     pause(0.1)
% end
% title('Segments')

%% Stokeslets Only
% Using exact velocities to recover forces using segment method, then
% calculate velocities everywhere using segment method (checking that
% segment method works).

% Segment points (left endpoints of segments)
y1seg = y1(segData(:,1)); % x-coords
y2seg = y2(segData(:,1)); % y-coords

% Using exact velocity from true solution at the segment points
[u1,u2,~] = permeablechannelexact(y1,y2,Da);
u = [u1, u2];
% Velocity in segment order
u1seg = u1(segData(:,1));
u2seg = u2(segData(:,1));

% Calculate forces using segment method
[fseg_temp] = st_segments_velocitytoforce([y1,y2],u,mu,segData);
f1seg = fseg_temp(:,1);
f2seg = fseg_temp(:,2);

% Calculate Channel Velocities (checking method)
[U_st_g] = st_segments_forcetovelocity(x,[y1seg,y2seg],fseg_temp,mu,segData);

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

% Step 1: Find g using u = utrue on solid bd and u = 0 on permeable bd
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

[gseg] = sd_segments_velocitytogforce([y1,y2],[u1b,u2b],normals,mu,segData);
g1seg = gseg(:,1);
g2seg = gseg(:,2);

% Step 2: Recover boundary velocity on permeable boundary using g
y1b = y1(idx); % boundary points of permeable section of boundary
y2b = y2(idx);

% Recover missing boundary velocity pieces on permeable section
[useg_perm] = sd_segments_gforcetovelocity([y1b,y2b],[y1,y2],gseg,normals,mu,segData);
useg_perm1 = useg_perm(:,1);
useg_perm2 = useg_perm(:,2);

u1b_seg = u1b;
u2b_seg = u2b;

% Replace the 0s with the newly computed velocities in permeable region
u1b_seg(idx) = useg_perm1;
u2b_seg(idx) = useg_perm2;

% Step 3: Use boundary velocities to solve for forces on boundary
ub_seg = [u1b_seg, u2b_seg];
[fseg] = st_segments_velocitytoforce([y1,y2],ub_seg,mu,segData);

% Step 4: Use forces on boundary to calculate velocities throughout channel
[useg_g] = st_segments_forcetovelocity(x,[y1seg,y2seg],fseg,mu,segData);

%% Figures

% Exact vs recovered boundary velocity
sk = 1;
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

% %% Check segments against MRS 
% 
% Ny1 = 45; % Number of y-points on top/bottom boundary 
% Ny2 = 15;  % Number of y-points on inlet/outlet
% 
% Nseg = 2*(Ny1+Ny2); % Total number of segments on the boundary 
% % Discretization step
% ds_x = (xmax - xmin)/(Ny1-1);
% %ds_y = (ymax - ymin)/(Ny2+1);
% ds_y = (ymax - ymin)/( ceil((ymax - ymin)/ds_x));
% 
% % Regularization parameter
% ep = ds_y;
% beta_value = Da*ep*C;
% 
% % Discretize Channel Walls  
% % Using midpoint to avoid the corners
% stb = xmin+ds_x/2:ds_x:xmax-ds_x/2; % top and bottom
% stb = stb';
% slr = ymin+ds_y/2:ds_y:ymax-ds_y/2; % left and right
% slr = slr';
% % Trapezoid
% % stb = xmin+ds_x:ds_x:xmax-ds_x; % top and bottom
% % stb = stb';
% % slr = ymin:ds_y:ymax; % left and right
% % slr = slr';
% 
% % Define coordinates on each wall (x-coord,y-coord)
% y1_top = stb;   y2_top = ymax*ones(size(stb)); % top wall coordinates (y1_top,y2_top)
% y1_bot = stb;   y2_bot = ymin*ones(size(stb)); % bottom wall coordinates (y1_bot,y2_bot)
% y1_left = xmin*ones(size(slr));   y2_left = slr; % left wall coordinates (y1_left,y2_left)
% y1_right = xmax*ones(size(slr));   y2_right = slr; % right wall coordinates (y1_right,y2_right)
% y1 = [y1_top; y1_bot; y1_left; y1_right]; %x-coordinates of all boundary points
% y2 = [y2_top; y2_bot; y2_left; y2_right]; %y-coordinates of all boundary points
% 
% % Trapezoid
% %wt = [ds_x/2; ds_x*ones(length(y1_top)-1,1); ds_y*ones(size(y1_bot)); ds_x*ones(size(y1_left)); ds_y*ones(length(y1_right)-1,1); ds_y/2];
% % Midpoint
% wt = [ds_x*ones(length(y1_top),1); ds_y*ones(size(y1_bot)); ds_x*ones(size(y1_left)); ds_y*ones(length(y1_right),1)];
% 
% % Indices of permeable region (where the boundary velocity is unknown)
% idx = find(-perm_factor*Lx<y1 & y1<perm_factor*Lx);
% solidI = find(y1 == xmin | y1 == xmax);
% beta = zeros(length(y1),1);
% beta(idx) = beta_value;
% 
% % Define normal vectors
% normals_top = zeros(length(stb),2); % unit normals for top:
% normals_top(:,2) = 1;
% normals_bot = zeros(length(stb),2); % unit normals for bottom
% normals_bot(:,2) = -1;
% normals_left = zeros(length(slr),2); % unit normals for left  wall 
% normals_left(:,1) = -1;
% normals_right = zeros(length(slr),2); % unit normals for right wall 
% normals_right(:,1) = 1;
% % normals on full boundary
% normals = [normals_top; normals_bot; normals_left; normals_right];
% % normals for top and bottom only
% normals_tb = [normals_top; normals_bot];
% 
% % Using exact velocity from true solution at the nodes
% [u1,u2,~] = permeablechannelexact(y1,y2,Da);
% 
% [fcheck] = RegStokeslets2D_velocitytoforce([y1 y2],[y1 y2],[u1 u2],ep,mu,1,wt);
% f1mrs = fcheck(:,1);
% f2mrs = fcheck(:,2);
% 
% u1b = u1; %x-coordinates of all boundary velocities
% u2b = u2; %y-coordinates of all boundary velocities
% u1b(idx) = 0; % set boundary velocities on permeable section to zero
% u2b(idx) = 0;
% 
% [gcheck] = RegStokeslets2D_velocityto_gforce_permeable([y1,y2],[y1,y2],...
%     [u1b,u2b], ep, mu, 1, idx, beta, normals, wt);
% g1mrs = gcheck(:,1);
% g2mrs = gcheck(:,2);
% 
% % Recover velocity on full boundary for comparison
% y1b = y1(idx); % boundary points of permeable section of boundary
% y2b = y2(idx);
% 
% % Recover missing velocity pieces on permeable section
% [ucheck_perm] = RegStokeslets2D_gtovelocity([y1,y2], gcheck, [y1b,y2b], ep, mu, 1, beta, normals, wt);
% ucheck_perm1 = ucheck_perm(:,1);
% ucheck_perm2 = ucheck_perm(:,2);
% 
% % Replace the placeholder zeros in u1,u2 with the velocities
% u1b(idx) = ucheck_perm1;
% u2b(idx) = ucheck_perm2;

%% Figures
% % Plot recovered vs exact boundary velocities
% sk = 1;
% figure
% subplot(2,1,1);
% hold on;
% plot(y1,y2,'k.')
% hold on
% quiver(y1(1:sk:end),y2(1:sk:end),u1b(1:sk:end)*scalefactor,u2b(1:sk:end)*scalefactor,'b','AutoScale','off')
% %quiver([y1_left; y1_right],[y2_left; y2_right],u1(203:end),u2(203:end),'b')
% title('Recovered Bd Velocities on Boundaries using MRS method')
% set(gca,'Fontsize',14)
% hold off;
% axis equal;
% axis([-Lx-1,Lx+1,-1.5,1.5])
% subplot(2,1,2);
% plot(y1,y2,'k.')
% hold on
% quiver(y1(1:sk:end),y2(1:sk:end),u1(1:sk:end)*scalefactor, u2(1:sk:end)*scalefactor,'r','Autoscale','off')
% hold off
% title('Exact Boundary Velocities')
% set(gca,'Fontsize',14)
% axis equal
% axis([-Lx-1,Lx+1,-1.5,1.5])

%%
% % Comparing Segments and MRS forces (Stokeslets)
% sk = 1;
% scaleFactor = 1.5;
% figure;
% subplot(2,1,1)
% hold on
% quiver(y1seg(1:sk:end),y2seg(1:sk:end),f1seg(1:sk:end),f2seg(1:sk:end),scaleFactor,'b')
% quiver(y1seg(1:sk:end),y2seg(1:sk:end),u1seg(1:sk:end),u2seg(1:sk:end),scaleFactor,'r')
% plot(y1seg(1:sk:end),y2seg(1:sk:end),'ro')
% axis equal
% title('Segments - Forces calculated from exact bd velocity')
% axis([-4,4,-1.5,1.5])
% subplot(2,1,2)
% hold on
% sk = 1;
% quiver(y1(1:sk:end),y2(1:sk:end),f1mrs(1:sk:end),f2mrs(1:sk:end),scaleFactor,'b')
% quiver(y1(1:sk:end),y2(1:sk:end),u1(1:sk:end),u2(1:sk:end),scaleFactor,'r')
% plot(y1(1:sk:end),y2(1:sk:end),'ro')
% axis equal
% title('MRS - Forces calculated from exact bd velocity')
% axis([-4,4,-1.5,1.5])
% 
% % Comparing Segment and MRS forces (source doublet)
% sk = 1;
% scaleFactor = 1.5;
% figure;
% subplot(2,1,1)
% hold on
% quiver(y1seg(1:sk:end),y2seg(1:sk:end),g1seg(1:sk:end),g2seg(1:sk:end),scaleFactor,'b')
% quiver(y1seg(1:sk:end),y2seg(1:sk:end),u1seg(1:sk:end),u2seg(1:sk:end),scaleFactor,'r')
% plot(y1seg(1:sk:end),y2seg(1:sk:end),'ro')
% axis equal
% title('Segments - SD forces calculated from exact bd velocity')
% axis([-4,4,-1.5,1.5])
% subplot(2,1,2)
% hold on
% sk = 1;
% quiver(y1(1:sk:end),y2(1:sk:end),g1mrs(1:sk:end),g2mrs(1:sk:end),scaleFactor,'b')
% quiver(y1(1:sk:end),y2(1:sk:end),u1(1:sk:end),u2(1:sk:end),scaleFactor,'r')
% plot(y1(1:sk:end),y2(1:sk:end),'ro')
% axis equal
% title('MRS - SD forces calculated from exact bd velocity')
% axis([-4,4,-1.5,1.5])

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
function [errors] = calculate_errors(u_error,v_error,xx1,xx2,dx,dy, method_name)

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
