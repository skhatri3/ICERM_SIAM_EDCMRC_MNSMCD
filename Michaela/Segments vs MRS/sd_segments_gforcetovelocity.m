function [uperm] = sd_segments_gforcetovelocity(x,y,g,normals,mu,segData)
%
% Calculates the velocity at observation points given the forces at the
% segments points using Stokeslet Segment method (2D). y = [y1 y2] and u
% should NOT be in the order of segData. This utilizes blob phi.

% Developed by Michaela Kubacki March 2026

% Inputs
%       x = (x1,x2) observation points
%       y = (y1,y2) source points (target points are source points)
%       g = (g1,g2) source doublet forces
%       normals = (outward) unit normal vectors
%       segData is a matrix with Nseg rows and at least 4 columns
%           segData(k,1) = index of the kth segment's start point
%           segData(k,2) = index of the kth segment's end point
%           segData(k,3) = epsilon value for the kth segment
%           segData(k,4) = beta value at segment k starting point
%           segData(k,5) = beta value at segment k ending point

% Unpack input
y1 = y(:,1);
y2 = y(:,2);

g1 = g(:,1);
g2 = g(:,2);

x1 = x(:,1);
x2 = x(:,2);

% Put in segment order
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

Nseg = length(y1_start);
Neval = length(x1);

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
    Norm1 = normalseg(k,1);
    Norm2 = normalseg(k,2);
    NdotXY = Norm1.*XY1 + Norm2.*XY2;
    
    % Calculate T and Q terms
    y = [y1_start(k) y2_start(k); y1_end(k) y2_end(k)];
    ep = segData(k,3);
    [T02,T04,T06,T12,T14,T16,T24,T26] = sd_segment_terms(x,y,ep);
    
    % Matrix blocks for u = Cfj + Dfk
    
    D11 = -((T12 + ep.^2.*T14).*Norm1.*Norm1 - NdotXY.*(4*ep.^2.*T16 + 2*T14).*XY1.*Norm1 - NdotXY.*(4*ep.^2.*T26 + 2*T24).*V1.*Norm1);
    D22 = -((T12 + ep.^2.*T14).*Norm2.*Norm2 - NdotXY.*(4*ep.^2.*T16 + 2*T14).*XY2.*Norm2 - NdotXY.*(4*ep.^2.*T26 + 2*T24).*V2.*Norm2);
    D12 = -((T12 + ep.^2.*T14).*Norm1.*Norm2 - NdotXY.*(4*ep.^2.*T16 + 2*T14).*XY1.*Norm2 - NdotXY.*(4*ep.^2.*T26 + 2*T24).*V1.*Norm2);
    D21 = -((T12 + ep.^2.*T14).*Norm1.*Norm2 - NdotXY.*(4*ep.^2.*T16 + 2*T14).*XY2.*Norm1 - NdotXY.*(4*ep.^2.*T26 + 2*T24).*V2.*Norm1);

    C11 = -((T02 + ep.^2.*T04).*Norm1.*Norm1 - NdotXY.*(4*ep.^2.*T06 + 2*T04).*XY1.*Norm1 - NdotXY.*(4*ep.^2.*T06 + 2*T04).*V1.*Norm1) - D11;
    C22 = -((T02 + ep.^2.*T04).*Norm2.*Norm2 - NdotXY.*(4*ep.^2.*T06 + 2*T04).*XY2.*Norm2 - NdotXY.*(4*ep.^2.*T06 + 2*T04).*V2.*Norm2) - D22;
    C12 = -((T02 + ep.^2.*T04).*Norm1.*Norm2 - NdotXY.*(4*ep.^2.*T06 + 2*T04).*XY1.*Norm2 - NdotXY.*(4*ep.^2.*T06 + 2*T04).*V1.*Norm2) - D12;
    C21 = -((T02 + ep.^2.*T04).*Norm1.*Norm2 - NdotXY.*(4*ep.^2.*T06 + 2*T04).*XY2.*Norm1 - NdotXY.*(4*ep.^2.*T06 + 2*T04).*V2.*Norm1) - D21;
    
    %Calculate velocity due to gforces on current segment, add to total
    %velocity
    u1(:) = u1(:) + L*(beta_start(k)*(C11*g1_start(k) + C12*g2_start(k)) + beta_end(k)*(D11*g1_end(k) + D12*g2_end(k)));
    u2(:) = u2(:) + L*(beta_start(k)*(C21*g1_start(k) + C22*g2_start(k)) + beta_end(k)*(D21*g1_end(k) + D22*g2_end(k)));

end

% Rescale
u1 = u1/(2*pi*mu);
u2 = u2/(2*pi*mu);

% Repack output
uperm = [u1 u2];