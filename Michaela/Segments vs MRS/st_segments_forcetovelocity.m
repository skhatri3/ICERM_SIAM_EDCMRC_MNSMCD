function [u] = st_segments_forcetovelocity(x,y,f,mu,segData,blob_num)
%
% Calculates the velocity at observation points given the forces at the
% segments points using Stokeslet Segment method (2D). y = [y1 y2] and f
% are not assumed to be in segment order. If they are in segment order,
% then segData should be segData(:,1) = (1:Nseg)'; segData(:,2) =
% [(2:Nseg)'; 1]. The boundary must be
% closed, with a cyclic segment ordering (so the end of the last segment is
% the starting point of the first segment).

% Developed by Michaela Kubacki and Ricardo Cortez March 2026

% Inputs
%       x = (x1,x2) observation points
%       y = (y1,y2) source points
%       f = (f1,f2) forces at the source points
%       segData is a matrix with Nseg rows and at least 3 columns
%           segData(k,1) = index of the kth segment's start point
%           segData(k,2) = index of the kth segment's end point
%           segData(k,3) = epsilon value for the kth segment
% Output: the velocity, u = [u1, u2] in segment order

% Unpack input
y1 = y(:,1);
y2 = y(:,2);

f1 = f(:,1);
f2 = f(:,2);

x1 = x(:,1);
x2 = x(:,2);

% Put in segment order
% Segment starting points
y1_start = y1(segData(:,1));
y2_start = y2(segData(:,1));
f1_start = f1(segData(:,1));
f2_start = f2(segData(:,1));

% Segment ending points
y1_end = y1(segData(:,2));
y2_end = y2(segData(:,2));
f1_end = f1(segData(:,2));
f2_end = f2(segData(:,2));

Nseg = length(y1_start); % Number of segments
Neval = length(x1); % Number of evaluation points

% Initialize variables
u1 = zeros(Neval,1);
u2 = zeros(Neval,1);
V1 = zeros(Neval,1);
V2 = zeros(Neval,1);

%% Loop over segments to calculate velocity at the evaluation points
for k = 1:Nseg
    % XY = xl - yk (differences between eval points and segment points)
    XY1 = x1(:) - y1_start(k); % x-coords of x-yk
    XY2 = x2(:) - y2_start(k); % y-coords of x-yk
    % V = yj - yk (tangent vector)
    V1(:) = (y1_start(k)-y1_end(k));
    V2(:) = (y2_start(k)-y2_end(k));
    L = sqrt(V1(1)^2+V2(1)^2); % length of current segment
    
    % Calculate T and Q terms
    y = [y1_start(k) y2_start(k); y1_end(k) y2_end(k)];
    ep = segData(k,3);    
    [H10, H11, H20, H21, H22, H23, ~, ~, ~, ~, ~] = seg_reg_fncs(x,y,ep,blob_num);

    % Stokeslet segment Matrix for ust = Afstart + Bfend

    B11 = (H11 + H21.*XY1.*XY1 + H22.*(XY1.*V1 + V1.*XY1) + H23.*V1.*V1);
    B22 = (H11 + H21.*XY2.*XY2 + H22.*(XY2.*V2 + V2.*XY2) + H23.*V2.*V2);
    B12 = (H21.*XY1.*XY2 + H22.*(XY1.*V2 + XY2.*V1)+ H23.*V1.*V2);

    A11 = (H10 + H20.*XY1.*XY1 + H21.*(XY1.*V1 +V1.*XY1) + H22.*V1.*V1) - B11;
    A22 = (H10 + H20.*XY2.*XY2 + H21.*(XY2.*V2 +V2.*XY2) + H22.*V2.*V2) - B22;
    A12 = (H20.*XY1.*XY2 + H21.*(XY1.*V2 +V1.*XY2) + H22.*V1.*V2) - B12;
    
    %Calculate velocity due to forces on current segment, add to total
    %velocity
    u1(:) = u1(:) + L*(A11*f1_start(k) + A12*f2_start(k) + B11*f1_end(k) + B12*f2_end(k));
    u2(:) = u2(:) + L*(A12*f1_start(k) + A22*f2_start(k) + B12*f1_end(k) + B22*f2_end(k));

end

% Rescale
u1 = u1/(mu);
u2 = u2/(mu);

% Repack output
u = [u1 u2];