function [uperm] = sd_segments_gforcetovelocity_beta_linear(x,y,gseg,normals,mu,segData,blob_num)
%
% Calculates the velocity at observation points given the forces at the
% segments points using Stokeslet Segment method (2D). y = [y1 y2] and u
% not assumed to be in the order of segData. If they are in segment order,
% then segData should be segData(:,1) = (1:Nseg)'; segData(:,2) =
% [(2:Nseg)'; 1]. The boundary must be closed, with a cyclic segment
% ordering (so the end of the last segment is the starting point of the
% first segment).

% Developed by Michaela Kubacki March 2026

% Inputs
%       x = (x1,x2) observation/evaluation points
%       y = (y1,y2) source points (target points are source points)
%       gseg = (g1,g2) source doublet forces (already in segment order)
%       normals = (outward) unit normal vectors
%       segData is a matrix with Nseg rows and at least 4 columns
%           segData(k,1) = index of the kth segment's start point
%           segData(k,2) = index of the kth segment's end point
%           segData(k,3) = epsilon value for the kth segment
%           segData(k,4) = beta value at segment k starting point
%           segData(k,5) = beta value at segment k ending point
%       blob_num = blob choice (1 for phi and 2 for psi)
%
% Output is the velocity uperm at the observation points x.

% Unpack input
y1  = y(:,1);
y2 = y(:,2);
g1 = gseg(:,1);
g2 = gseg(:,2);

x1 = x(:,1);
x2 = x(:,2);

% Specify variables in terms of segments

% Normals on segments
normalseg = normals(segData(:,1),:);

% Segment starting points
y1_start = y1(segData(:,1));
y2_start = y2(segData(:,1));
g1_start = g1(segData(:,1));
g2_start = g2(segData(:,1));
beta_start = segData(:,4);

% Segment ending points
y1_end = y1(segData(:,2));
y2_end = y2(segData(:,2));
g1_end = g1(segData(:,2));
g2_end = g2(segData(:,2));
beta_end = segData(:,5);

Nseg = length(y1_start); % Number of segments
Neval = length(x1); % Number of evaluation points

% Beta is linear on each segment
beta = beta_start;
betaHat = beta_end - beta_start;

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
    % V = yj - yk (tangent vector, start - finish)
    V1(:) = (y1_start(k)-y1_end(k));
    V2(:) = (y2_start(k)-y2_end(k));
    L = sqrt(V1(1)^2+V2(1)^2); % length of current segment
    Norm1 = normalseg(k,1); % x-coords of current segment normal
    Norm2 = normalseg(k,2); % y-coords of current segment normal
    NdotXY = Norm1.*XY1 + Norm2.*XY2;
    
    % Calculate segment functions for Source Doublet
    y = [y1_start(k) y2_start(k); y1_end(k) y2_end(k)];
    ep = segData(k,3);    
    [~, ~, ~, ~, ~, ~, S10, S11, S20, S21, S22, S12, S23] = seg_reg_fncs(x,y,ep,blob_num);
    
    % Matrix Blocks for segment solution
    D11 = beta(k)*(Norm1.*S11.*Norm1 + NdotXY.*S21.*XY1.*Norm1 + NdotXY.*S22.*V1.*Norm1) + ...
            betaHat(k)*(Norm1.*S12.*Norm1 + NdotXY.*S22.*XY1.*Norm1 + NdotXY.*S23.*V1.*Norm1);
    D22 = beta(k)*(Norm2.*S11.*Norm2 + NdotXY.*S22.*XY2.*Norm2 + NdotXY.*S23.*V2.*Norm2) + ...
            betaHat(k)*(Norm2.*S12.*Norm2 + NdotXY.*S21.*XY2.*Norm2 + NdotXY.*S22.*V2.*Norm2);
    D12 = beta(k)*(Norm1.*S11.*Norm2 + NdotXY.*S21.*XY1.*Norm2 + NdotXY.*S22.*V1.*Norm2) + ...
            betaHat(k)*(Norm1.*S12.*Norm2 + NdotXY.*S22.*XY1.*Norm2 + NdotXY.*S23.*V1.*Norm2);
    D21 = beta(k)*(Norm1.*S11.*Norm2 + NdotXY.*S21.*XY2.*Norm1 + NdotXY.*S22.*V2.*Norm1) + ...
            betaHat(k)*(Norm1.*S12.*Norm2 + NdotXY.*S22.*XY2.*Norm1 + NdotXY.*S23.*V2.*Norm1);

    C11 = beta(k)*(Norm1.*S10.*Norm1 + NdotXY.*S20.*XY1.*Norm1 + NdotXY.*S21.*V1.*Norm1) + ...
            betaHat(k)*(Norm1.*S11.*Norm1 + NdotXY.*S21.*XY1.*Norm1 + NdotXY.*S22.*V1.*Norm1)- D11;
    C22 = beta(k)*(Norm2.*S10.*Norm2 + NdotXY.*S20.*XY2.*Norm2 + NdotXY.*S21.*V2.*Norm2) + ...
            betaHat(k)*(Norm2.*S11.*Norm2 + NdotXY.*S21.*XY2.*Norm2 + NdotXY.*S22.*V2.*Norm2)- D22;
    C12 = beta(k)*(Norm1.*S10.*Norm2 + NdotXY.*S20.*XY1.*Norm2 + NdotXY.*S21.*V1.*Norm2) + ...
            betaHat(k)*(Norm1.*S11.*Norm2 + NdotXY.*S21.*XY1.*Norm2 + NdotXY.*S22.*V1.*Norm2)- D12;
    C21 = beta(k)*(Norm1.*S10.*Norm2 + NdotXY.*S20.*XY2.*Norm1 + NdotXY.*S21.*V2.*Norm1) + ...
            betaHat(k)*(Norm1.*S11.*Norm2 + NdotXY.*S21.*XY2.*Norm1 + NdotXY.*S22.*V2.*Norm1)- D21;
    
    %Calculate velocity due to gforces on current segment, add to total
    %velocity
    u1(:) = u1(:) - L*(C11*g1_start(k) + C12*g2_start(k) + D11*g1_end(k) + D12*g2_end(k));
    u2(:) = u2(:) - L*(C21*g1_start(k) + C22*g2_start(k) + D21*g1_end(k) + D22*g2_end(k));
end

% Rescale
u1 = u1/(mu);
u2 = u2/(mu);

% Repack output
uperm = [u1 u2];