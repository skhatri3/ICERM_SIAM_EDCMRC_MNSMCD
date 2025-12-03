function [g] = RegStokeslets2D_velocityto_gforce_permeable(y,...
    x,u,ep,mu,blob_num, I, beta, normal, wt)


% Computes intermediate "g" forces for permeable membrane 
% using the Method of Regularized Stokeslets. 
% Based on Cortez, SIAM J. Sci Comput. 2001 and Cortez, Fluids 2021 

% Developed by Shilpa Khatri, Ricardo Cortez, Brittany Leathers, and
% Michaela Kubacki July 2024.

% This function requires access to function reg_funcs_withdoublet.m

% Modified December 2025 to incorporate a quad weight vector.

% Inputs:
%       y = (y1,y2) source points 
%       x = (x1,x2) target points 
%       u = (u1,u2) velocity at those target points 
%       ep is the width of the regularization 
%       mu is the viscosity 
%       blob_num = blob choice number
%       I = Indices of x (target) rows corresponding to the permeable region
%       beta = numerical permeability coefficient for all source points
%       beta(i) should be 0 if it is not a permeable section
%       normal = unit normal vectors for all source points
%       wt = quadrature weights for all source points 

% ACCURACTY NOTE: If inf norm of condition number greater than 10^r (r>0)
%                then losing r digits of accuracy in matrix solve. Warning
%                message will display if r >= 8.

N = size(y,1); % number of source points 
M = size(x,1); % number of target points 

% unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2);
u1 = u(:,1);
u2 = u(:,2); 
x1 = x(:,1); 
x2 = x(:,2); 

% building matrix 

Beta = zeros(M,N);
Norm1 = zeros(M,N);
Norm2 = zeros(M,N);
XY1 = zeros(M,N);
XY2 = zeros(M,N);

%loop over source points    
for k = 1:N 

    % distance between target and source points 
    XY1(:,k) = x1(:) - y1(k); 
    XY2(:,k) = x2(:) - y2(k);  

    % Beta matrix: column k has beta_k
    Beta(:, k)=beta(k);

    % Normal matrices needed: 
    Norm1(:,k)=normal(k,1); %kth column is all nk1
    Norm2(:,k)=normal(k,2); %kth column is all nk2    
   
end


% Stokeslet and Source Doublet solution pieces for blob choice
R2 = XY1.^2 + XY2.^2 + ep^2; 
R = sqrt( R2 ); 

[H1, H2, S1, S2, ~] = reg_fncs_withdoublet(ep, R, blob_num); 

% Stokeslet Matrix
M11 = (H1 + H2.*XY1.*XY1).*wt; 
M22 = (H1 + H2.*XY2.*XY2).*wt; 
M12 = (H2.*XY1.*XY2).*wt; % Note that M21 = M12

% Assemble Stokeslet matrix and rescale
Mat = [M11 M12; M12 M22]/(mu); 

% Doublet Matrix
NormXY = Norm1.*XY1 + Norm2.*XY2;
D11 = -Beta.*Norm1 .* (S1.*Norm1+S2.*NormXY.*XY1) .*wt;
D12 = -Beta.*Norm2 .* (S1.*Norm1+S2.*NormXY.*XY1) .*wt;
D21 = -Beta.*Norm1 .* (S1.*Norm2+S2.*NormXY.*XY2) .*wt;
D22 = -Beta.*Norm2 .* (S1.*Norm2+S2.*NormXY.*XY2) .*wt;

% Assemble Double matrix and rescale
D = [D11 D12; D21 D22]/(mu);

% Zero out rows corresponding target points on permeable membrane
D(I,:)=zeros(length(I),length(D(1,:))); % corresponding to u1 components
D(I+M,:)=zeros(length(I), length(D(1,:))); % corresponding to u2 components

% Solving for force when given a velocity 
uu = [u1;u2];

% Conditioning check
disp(['log10 of condition number of Mat+D: ', num2str(log10(condest(Mat+D)))])
disp(['condition number of Mat+D: ', num2str(condest(Mat+D))])
disp(['size of Mat+D: ', num2str(size(Mat+D))])
log10cond = log10(condest(Mat+D));
if log10cond >= 8
    disp(['Warning: High condition number. Losing roughly ', num2str(log10cond), ' digits of accuracy.'])
end

gg = (Mat + D) \uu; 
g1 = gg(1:N); 
g2 = gg(N+1:end); 

% repacking output 
g = [g1, g2]; 
