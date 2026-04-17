function [g] = sd_segments_velocitytogforce_beta_linear(y,ub,normals,mu,segData,blob_num)
%
%
% Calculates intermediate "g" forces at the segment points on a permeable
% membrane boundary given the velocities at the segment points on the solid
% boundary using Source Doublet Segment method (2D). y = [y1 y2] and u
% should NOT be in the order of segData. If they are in segment order, then
% segData should be segData(:,1) = (1:Nseg)'; segData(:,2) = [(2:Nseg)';
% 1]. The boundary must be closed, with a cyclic segment ordering (so the
% end of the last segment is the starting point of the first segment).

% Developed by Michaela Kubacki March 2026
%
% Inputs
%       y = (y1,y2) source points (target points are source points)
%       ub = (ub1,ub2) velocities at source points, set to zero on the
%            permeable portion of the boundary
%       normals = (outward) unit normal vectors
%       segData is a matrix with Nseg rows and at least 4 columns
%           segData(k,1) = index of the kth segment's start point
%           segData(k,2) = index of the kth segment's end point
%           segData(k,3) = epsilon value for the kth segment
%           segData(k,4) = beta value at segment k starting point
%           segData(k,5) = beta value at segment k ending point
%       blob_num = blob choice (1 for phi and 2 for psi)
%
% Output is source doublet distribution, g, at source points, y (in segment
% order).
%
% Unpack input
y1 = y(:,1);
y2 = y(:,2);

% Put in segment order
y1seg = y1(segData(:,1));
y2seg = y2(segData(:,1));
beta_start = segData(:,4);

y1seg_end = y1(segData(:,2));
y2seg_end = y2(segData(:,2));
beta_end = segData(:,5);

ub1 = ub(:,1);
ub2 = ub(:,2);
ub1seg = ub1(segData(:,1));
ub2seg = ub2(segData(:,1));

normalseg = normals(segData(:,1),:);


% Find indeces of target points (segment start points) where beta is
% nonzero. These are points in the permeable region
idx = find(segData(:,4) ~=0);

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
epMat = XY1; % epsilon values of segments
Norm1 = XY1; Norm2 = XY1; % x and y coords of segment normals
beta = XY1; % beta values of segments
betaHat = XY1;

% Initializing segment reg funcs matrices
H10 = XY1; H11 = XY1; H20 = XY1; H21 = XY1; H22 = XY1; H23 = XY2;
S10 = XY1; S11 = XY1; S12 = XY1; S20 = XY1; S21 = XY1; S22 = XY1; S23 = XY1;


%% Loop over segments to create Nseg x Nseg block matrices
% Mat(l,k) = value of quantity you get using eval point xl on segment k
for k = 1:Nseg
    % XY = xl - yk (differences between eval points and segment points)
    XY1(:,k) = x1(:) - y1seg(k); % x-coords of x-yk
    XY2(:,k) = x2(:) - y2seg(k); % y-coords of x-yk
    % y = [left endpoint; right endpoint]
    y = [y1seg(k) y2seg(k); y1seg_end(k) y2seg_end(k)]; % start-end points of segment k
    V1Mat(:,k) = y(1,1)-y(2,1); % y1(k)-y1(k+1) x-coords of tangent vector
    V2Mat(:,k) = y(1,2)-y(2,2); % y2(k)-y2(k+1) y-coords of tangent vector
    LMat(:,k) = sqrt(V1Mat(:,k).^2 + V2Mat(:,k).^2); % length of segment k
    epMat(:,k) = segData(k,3); % epsilon value on segment k
    beta(:,k) = beta_start(k); % starting beta value of segment k
    betaHat(:,k) = beta_end(k)-beta_start(k); % difference in beta values
    Norm1(:,k) = normalseg(k,1); % x-coords of normal on segment k
    Norm2(:,k) = normalseg(k,2); % y-coords of normal on segment k
    % Calculate segment reg funcs needed for segment matrix formula
    [H10(:,k), H11(:,k), H20(:,k), H21(:,k), H22(:,k), H23(:,k), S10(:,k), ...
        S11(:,k), S20(:,k), S21(:,k), S22(:,k),S12(:,k),S23(:,k)] = seg_reg_fncs(x,y,epMat(:,k),blob_num);
end


%% Build matrix
% Rows correspond to eval points, columns to segments.
NdotXY = Norm1.*XY1 + Norm2.*XY2;

% Stokeslet segment Matrix for ust = Afstart + Bfend

B11tmp = (H11 + H21.*XY1.*XY1 + H22.*(XY1.*V1Mat + V1Mat.*XY1) + H23.*V1Mat.*V1Mat).*LMat;
B22tmp = (H11 + H21.*XY2.*XY2 + H22.*(XY2.*V2Mat +V2Mat.*XY2) + H23.*V2Mat.*V2Mat).*LMat;
B12tmp = (H21.*XY1.*XY2 + H22.*(XY1.*V2Mat + XY2.*V1Mat)+ H23.*V1Mat.*V2Mat).*LMat; % B12 = B21

A11 = (H10 + H20.*XY1.*XY1 + H21.*(XY1.*V1Mat +V1Mat.*XY1) + H22.*V1Mat.*V1Mat).*LMat - B11tmp;
A22 = (H10 + H20.*XY2.*XY2 + H21.*(XY2.*V2Mat +V2Mat.*XY2) + H22.*V2Mat.*V2Mat).*LMat - B22tmp;
A12 = (H20.*XY1.*XY2 + H21.*(XY1.*V2Mat +V1Mat.*XY2) + H22.*V1Mat.*V2Mat).*LMat - B12tmp; % A12 = A21

% For closed curve, where first starting segment point = last segment
% ending point
B11 = circshift(B11tmp,1,2);
B22 = circshift(B22tmp,1,2);
B12 = circshift(B12tmp,1,2);

% Assemble Stokeslet Segment Matrix
Mst = [A11+B11 A12+B12; A12+B12 A22+B22]/mu;

% Source Doublet Matrix for usd = Cfstart + Dfend;, where f = beta*g

D11tmp = beta.*(S11.*Norm1.*Norm1 + NdotXY.*S21.*XY1.*Norm1 + NdotXY.*S22.*V1Mat.*Norm1).*LMat + ...
            betaHat.*(S12.*Norm1.*Norm1 + NdotXY.*S22.*XY1.*Norm1 + NdotXY.*S23.*V1Mat.*Norm1).*LMat;

D22tmp = beta.*(S11.*Norm2.*Norm2 + NdotXY.*S21.*XY2.*Norm2 + NdotXY.*S22.*V2Mat.*Norm2).*LMat + ...
            betaHat.*(S12.*Norm2.*Norm2 + NdotXY.*S22.*XY2.*Norm2 + NdotXY.*S23.*V2Mat.*Norm2).*LMat;

D12tmp = beta.*(S11.*Norm1.*Norm2 + NdotXY.*S21.*XY1.*Norm2 + NdotXY.*S22.*V1Mat.*Norm2).*LMat + ...
            betaHat.*(S12.*Norm1.*Norm2 + NdotXY.*S22.*XY1.*Norm2 + NdotXY.*S23.*V1Mat.*Norm2).*LMat;

D21tmp = beta.*(S11.*Norm2.*Norm1 + NdotXY.*S21.*XY2.*Norm1 + NdotXY.*S22.*V2Mat.*Norm1).*LMat + ...
            betaHat.*(S12.*Norm2.*Norm1 + NdotXY.*S22.*XY2.*Norm1 + NdotXY.*S23.*V2Mat.*Norm1).*LMat;


C11 = beta.*((S10.*Norm1.*Norm1 + NdotXY.*S20.*XY1.*Norm1 + NdotXY.*S21.*V1Mat.*Norm1).*LMat) + ...
        betaHat.*((S11.*Norm1.*Norm1 + NdotXY.*S21.*XY1.*Norm1 + NdotXY.*S22.*V1Mat.*Norm1).*LMat) - D11tmp;

C22 = beta.*((S10.*Norm2.*Norm2 + NdotXY.*S20.*XY2.*Norm2 + NdotXY.*S21.*V2Mat.*Norm2).*LMat) + ...
        betaHat.*((S11.*Norm2.*Norm2 + NdotXY.*S21.*XY2.*Norm2 + NdotXY.*S22.*V2Mat.*Norm2).*LMat)- D22tmp;

C12 = beta.*((S10.*Norm1.*Norm2 + NdotXY.*S20.*XY1.*Norm2 + NdotXY.*S21.*V1Mat.*Norm2).*LMat) + ... 
        betaHat.*((S11.*Norm1.*Norm2 + NdotXY.*S21.*XY1.*Norm2 + NdotXY.*S22.*V1Mat.*Norm2).*LMat)- D12tmp;

C21 = beta.*((S10.*Norm2.*Norm1 + NdotXY.*S20.*XY2.*Norm1 + NdotXY.*S21.*V2Mat.*Norm1).*LMat)  + ...
        betaHat.*((S11.*Norm2.*Norm1 + NdotXY.*S21.*XY2.*Norm1 + NdotXY.*S22.*V2Mat.*Norm1).*LMat)- D21tmp;


% For closed curve, where first starting segment point = last segment
% ending point
D11 = circshift(D11tmp,1,2);
D22 = circshift(D22tmp,1,2);
D12 = circshift(D12tmp,1,2);
D21 = circshift(D21tmp,1,2);

% Assemble Source Double Segment matrix
Msd = -[C11+D11 C12+D12; C21+D21 C22+D22]/mu;

% Zero out rows corresponding to eval points on permeable membrane to preserve
% Mst*g = 0 on permeable membrane (known boundary velocity is the sum of
% contributions from St and SD, unknown boundary velocity is entirely from
% SD). 

Msd(idx,:) = zeros(length(idx),length(Msd(1,:))); % corresponding to u1 components
Msd(idx + Nseg,:) = zeros(length(idx),length(Msd(1,:))); % corresponding to u2 components

%% Solve for gforce given boundary velocity
Useg = [ub1seg; ub2seg];

G = (Mst + Msd) \ Useg;
g1 = G(1:Nseg);
g2 = G(Nseg+1:end);

% Repacking output
g = [g1 g2];

