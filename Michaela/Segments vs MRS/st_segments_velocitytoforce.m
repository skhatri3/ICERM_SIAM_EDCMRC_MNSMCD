function [f] = st_segments_velocitytoforce(y,u,mu,segData,blob_num)
%
% Calculates the forces at the segment points given the velocities at the
% segments points using Stokeslet Segment method (2D). y = [y1 y2] and u
% should NOT be in the order of segData. y = [y1 y2] and u are not assumed
% to be in segment order. If they are in segment order, then segData should
% be segData(:,1) = (1:Nseg)'; segData(:,2) = [(2:Nseg)'; 1]. The boundary
% must be closed, with a cyclic segment ordering (so the end of the last
% segment is the starting point of the first segment).

% Developed by Michaela Kubacki and Ricardo Cortez March 2026
%
% Inputs
%       y = (y1,y2) source points (in this case, source = target points)
%       u = (u1,u2) velocities at source points
%       segData is a matrix with Nseg rows and at least 3 columns
%           segData(k,1) = index of the kth segment's start point
%           segData(k,2) = index of the kth segment's end point
%           segData(k,3) = epsilon value for the kth segment
%       blob_num = blob choice (1 for phi and 2 for psi)
%
% Output force f, evaluated at source points, y (in segment order).
%
% Unpack input
y1 = y(:,1);
y2 = y(:,2);

u1 = u(:,1);
u2 = u(:,2);

% Put in segment order
y1seg = y1(segData(:,1));
y2seg = y2(segData(:,1));

y1seg_end = y1(segData(:,2));
y2seg_end = y2(segData(:,2));

u1seg = u1(segData(:,1));
u2seg = u2(segData(:,1));

% Eval points are segment endpoints
x1 = y1seg;
x2 = y2seg;
x = [x1, x2];

Nseg = length(y1seg);

% Initializing matrices of size Nseg x Nseg
XY1 = zeros(Nseg,Nseg); % first component of x-yk
XY2 = XY1; % second component of x-yk
V1Mat = XY1; % first component of yl-y(k+1)
V2Mat = XY1; % second component of yl-y(k+1)
LMat = XY1; % lengths of segments
epMat = XY1; % epsilon values for segments
% Initialize matrices for segment reg funcs
H10 = XY1; H11 = XY1; H20 = XY1; H21 = XY1; H22 = XY1; H23 = XY2;


%% Loop over segments to create Nseg x Nseg block matrices
% Mat(l,k) = value of quantity you get using eval point xl on segment k
for k = 1:Nseg
    % XY = xl - yk (differences between eval points and segment points)
    XY1(:,k) = x1(:) - y1seg(k); % x-coords of x-yk
    XY2(:,k) = x2(:) - y2seg(k); % y-coords of x-yk
    % y = [left endpoint; right endpoint]
    y = [y1seg(k) y2seg(k); y1seg_end(k) y2seg_end(k)];
    V1Mat(:,k) = y(1,1)-y(2,1); % y1(k)-y1(k+1) x-coords of tangent vector
    V2Mat(:,k) = y(1,2)-y(2,2); % y2(k)-y2(k+1) y-coords of tangent vector
    LMat(:,k) = sqrt(V1Mat(:,k).^2 + V2Mat(:,k).^2);
    % Calculate segment reg functions needed for segment matrix formula
    epMat(:,k) = segData(k,3)*ones(Nseg,1);
    [H10(:,k), H11(:,k), H20(:,k), H21(:,k), H22(:,k), H23(:,k), ~, ~, ~, ~, ~] = seg_reg_fncs(x,y,epMat(:,k),blob_num);
end

%% Build matrix
% Rows correspond to eval points, columns to segments.

% Stokeslet segment Matrix for ust = Afstart + Bfend

B11tmp = (H11 + H21.*XY1.*XY1 + H22.*(XY1.*V1Mat + V1Mat.*XY1) + H23.*V1Mat.*V1Mat).*LMat;
B22tmp = (H11 + H21.*XY2.*XY2 + H22.*(XY2.*V2Mat +V2Mat.*XY2) + H23.*V2Mat.*V2Mat).*LMat;
B12tmp = (H21.*XY1.*XY2 + H22.*(XY1.*V2Mat + XY2.*V1Mat)+ H23.*V1Mat.*V2Mat).*LMat;

A11 = (H10 + H20.*XY1.*XY1 + H21.*(XY1.*V1Mat +V1Mat.*XY1) + H22.*V1Mat.*V1Mat).*LMat - B11tmp;
A22 = (H10 + H20.*XY2.*XY2 + H21.*(XY2.*V2Mat +V2Mat.*XY2) + H22.*V2Mat.*V2Mat).*LMat - B22tmp;
A12 = (H20.*XY1.*XY2 + H21.*(XY1.*V2Mat +V1Mat.*XY2) + H22.*V1Mat.*V2Mat).*LMat - B12tmp;

B11 = circshift(B11tmp,1,2);
B22 = circshift(B22tmp,1,2);
B12 = circshift(B12tmp,1,2);

M = [A11+B11 A12+B12; A12+B12 A22+B22]/(mu);

%% Solve linear system
U = [u1seg;u2seg];
F = M\U;
f1 = F(1:Nseg);
f2 = F(Nseg+1:end);
f = [f1 f2];
end