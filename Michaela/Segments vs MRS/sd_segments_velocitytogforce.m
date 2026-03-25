function [g] = sd_segments_velocitytogforce(y,ub,normals,mu,segData)
%
% Calculates intermediate "g" forces at the segment points on a permeable
% membrane boundary given the velocities at the segments points on the
% solid boundary using Source Doublet Segment method (2D). y = [y1 y2] and
% u should NOT be in the order of segData. The segments should be ordered
% in a way that each node is the endpoint of at most one segment. The
% boundary must be closed.

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
%
% Unpack input
y1 = y(:,1);
y2 = y(:,2);

% Put in segment order
y1seg = y1(segData(:,1));
y2seg = y2(segData(:,1));

y1seg_end = y1(segData(:,2));
y2seg_end = y2(segData(:,2));

ub1 = ub(:,1);
ub2 = ub(:,2);
ub1seg = ub1(segData(:,1));
ub2seg = ub2(segData(:,1));

normalseg = normals(segData(:,1),:);

% Find indeces of permeable segments (where beta is nonzero)
idx_seg_start = find(segData(:,4) ~= 0); 
idx_seg_end = find(segData(:,5) ~= 0);
idx_seg = union(idx_seg_start,idx_seg_end);

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
T02Mat = XY1; T12Mat = XY1; T22Mat = XY1; T32Mat = XY1;
T04Mat = XY1; T06Mat = XY1; 
T14Mat = XY1; T16Mat = XY1; 
T24Mat = XY1; T26Mat = XY1;
Q0Mat = XY1; Q1Mat = XY1;
epMat = XY1; beta_start = XY1; beta_end = XY1;
Norm1 = XY1; Norm2 = XY1;

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
    epMat(:,k) = segData(k,3)*ones(Nseg,1);
    beta_start(:,k) = segData(k,4)*ones(Nseg,1);
    beta_end(:,k) = segData(k,5)*ones(Nseg,1);
    Norm1(:,k) = normalseg(k,1);
    Norm2(:,k) = normalseg(k,2);
    % Calculate values needed for segment matrix formula
    [Q0,Q1,T02,T12,T22,T32] = st_segment_terms(x,y,segData(k,3));
    T02Mat(:,k) = T02;
    T12Mat(:,k) = T12;
    T22Mat(:,k) = T22;
    T32Mat(:,k) = T32;
    Q1Mat(:,k) = Q1;
    Q0Mat(:,k) = Q0;
    [~,T04,T06,~,T14,T16,T24,T26] = sd_segment_terms(x,y,segData(k,3));
    T04Mat(:,k) = T04;
    T06Mat(:,k) = T06;
    T14Mat(:,k) = T14;
    T16Mat(:,k) = T16;
    T24Mat(:,k) = T24;
    T26Mat(:,k) = T26;
end

NdotXY = Norm1.*XY1 + Norm2.*XY2;

%% Build matrix
% Rows correspond to eval points, columns to segments.

B11tmp = (epMat.^2.*T12Mat - Q1Mat + T12Mat.*XY1.*XY1 + T22Mat.*(XY1.*V1Mat + V1Mat.*XY1) + T32Mat.*V1Mat.*V1Mat).*LMat;
B22tmp = (epMat.^2.*T12Mat - Q1Mat + T12Mat.*XY2.*XY2 + T22Mat.*(XY2.*V2Mat +V2Mat.*XY2) + T32Mat.*V2Mat.*V2Mat).*LMat;
B12tmp = (T12Mat.*XY1.*XY2 + T22Mat.*(XY1.*V2Mat + XY2.*V1Mat)+ T32Mat.*V1Mat.*V2Mat).*LMat;

A11 = (epMat.^2.*T02Mat - Q0Mat + T02Mat.*XY1.*XY1 + T12Mat.*(XY1.*V1Mat +V1Mat.*XY1) + T22Mat.*V1Mat.*V1Mat).*LMat - B11tmp;
A22 = (epMat.^2.*T02Mat - Q0Mat + T02Mat.*XY2.*XY2 + T12Mat.*(XY2.*V2Mat +V2Mat.*XY2) + T22Mat.*V2Mat.*V2Mat).*LMat - B22tmp;
A12 = (T02Mat.*XY1.*XY2 + T12Mat.*(XY1.*V2Mat +V1Mat.*XY2) + T22Mat.*V1Mat.*V2Mat).*LMat - B12tmp;

% For closed curve, where first starting segment point = last segment
% ending point
B11 = circshift(B11tmp,1,2);
B22 = circshift(B22tmp,1,2);
B12 = circshift(B12tmp,1,2);

% Assemble Stokeslet Segment Matrix
Mst = [A11+B11 A12+B12; A12+B12 A22+B22]/(4*pi*mu);

% Note about beta: on segment k: u_sd = -L(C_k*(beta_start*f_start)+ D_k
% *(beta_end f_end))

D11tmp = ((T12Mat + epMat.^2.*T14).*Norm1.*Norm1 - NdotXY.*(4*epMat.^2.*T16Mat + 2*T14Mat).*XY1.*Norm1 - NdotXY.*(4*epMat.^2.*T26Mat + 2*T24Mat).*V1Mat.*Norm1).*LMat;
D22tmp = ((T12Mat + epMat.^2.*T14).*Norm2.*Norm2 - NdotXY.*(4*epMat.^2.*T16Mat + 2*T14Mat).*XY2.*Norm2 - NdotXY.*(4*epMat.^2.*T26Mat + 2*T24Mat).*V2Mat.*Norm2).*LMat;
D12tmp = ((T12Mat + epMat.^2.*T14).*Norm1.*Norm2 - NdotXY.*(4*epMat.^2.*T16Mat + 2*T14Mat).*XY1.*Norm2 - NdotXY.*(4*epMat.^2.*T26Mat + 2*T24Mat).*V1Mat.*Norm2).*LMat;
D21tmp = ((T12Mat + epMat.^2.*T14).*Norm1.*Norm2 - NdotXY.*(4*epMat.^2.*T16Mat + 2*T14Mat).*XY2.*Norm1 - NdotXY.*(4*epMat.^2.*T26Mat + 2*T24Mat).*V2Mat.*Norm1).*LMat;

C11 = beta_start.*(((T02Mat + epMat.^2.*T04Mat).*Norm1.*Norm1 - NdotXY.*(4*epMat.^2.*T06Mat + 2*T04Mat).*XY1.*Norm1 - NdotXY.*(4*epMat.^2.*T06Mat + 2*T04Mat).*V1Mat.*Norm1).*LMat - D11tmp);
C22 = beta_start.*(((T02Mat + epMat.^2.*T04Mat).*Norm2.*Norm2 - NdotXY.*(4*epMat.^2.*T06Mat + 2*T04Mat).*XY2.*Norm2 - NdotXY.*(4*epMat.^2.*T06Mat + 2*T04Mat).*V2Mat.*Norm2).*LMat - D22tmp);
C12 = beta_start.*(((T02Mat + epMat.^2.*T04Mat).*Norm1.*Norm2 - NdotXY.*(4*epMat.^2.*T06Mat + 2*T04Mat).*XY1.*Norm2 - NdotXY.*(4*epMat.^2.*T06Mat + 2*T04Mat).*V1Mat.*Norm2).*LMat - D12tmp);
C21 = beta_start.*(((T02Mat + epMat.^2.*T04Mat).*Norm1.*Norm2 - NdotXY.*(4*epMat.^2.*T06Mat + 2*T04Mat).*XY2.*Norm1 - NdotXY.*(4*epMat.^2.*T06Mat + 2*T04Mat).*V2Mat.*Norm1).*LMat - D21tmp);

D11tmp = beta_end.*D11tmp;
D22tmp = beta_end.*D22tmp;
D12tmp = beta_end.*D12tmp;
D21tmp = beta_end.*D11tmp;

% For closed curve, where first starting segment point = last segment
% ending point
D11 = circshift(D11tmp,1,2);
D22 = circshift(D22tmp,1,2);
D12 = circshift(D12tmp,1,2);
D21 = circshift(D21tmp,1,2);

Msd = -[C11+D11 C12+D12; C21+D21 C22+D22]/(2*pi*mu);

% Zero out rows corresponding to points on permeable membrane to preserve
% Mst*g = 0 on permeable membrane
% NOTE: Not sure if this is correct! Took the union of the left and right indeces
% where beta is nonzero. Otherwise plots were not symmetric!
Msd(idx_seg,:) = zeros(length(idx_seg),length(Msd(1,:))); 
Msd(idx_seg + Nseg,:) = zeros(length(idx_seg),length(Msd(1,:)));

%% Solve for gforce given boundary velocity
Useg = [ub1seg; ub2seg];

G = (Mst + Msd) \ Useg;
g1 = G(1:Nseg);
g2 = G(Nseg+1:end);

% Repacking output
g = [g1 g2];

