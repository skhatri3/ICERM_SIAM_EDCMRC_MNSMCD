function [u_beta] = RegStokeslets2D_gtovelocity_2eps(y,g,...
    x,ep1, ep2,mu,blob_num, beta, normal)

% Computes velocities on permeable membrane from intermediate force g
% using the Method of Regularized Stokeslets 
% Based on Cortez, SIAM J. Sci Comput. 2001 and Cortez, Fluids 2021 

% Developed by Shilpa Khatri, Ricardo Cortez, Brittany Leathers, and
% Michaela Kubacki
% July 2024 

%y = (y1,y2) source points
%g = (g1,g2) forces at those source points (not force density) Nx2
%%%%%x = (x1,x2) target points - these should be just the ones in the
%%%%%permeable part of the membrane
%u = (u1,u2) velocity evaluated at those target points 
%mu is the viscosity 
%ep is the width of the regularization (ep1 for stokeslet part, ep2 for
%doublet part)
%I: Indices of x (target) rows that correspond to Gamma_beta, the permeable region
%Beta gives the permeability coefficient for all source points., Nx1
    %beta(i) should be 0 if it is not a permeable section
%n is Nx2, giving unit normals for source points
%I2 is indicies of non-permeable parts of y1, but it has been taken out of
%the code

N = size(y,1); %number of source points 
M = size(x,1); %number of target points 

%unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2);
g1 = g(:,1); 
g2 = g(:,2); 
x1 = x(:,1); 
x2 = x(:,2); 


%initializing the velocity 
u1 = zeros(M,1);
u2 = zeros(M,1); 


%loop over source points    
for k = 1:N 

    %distance between target and source points 
    XY1 = x1(:) - y1(k); 
    XY2= x2(:) - y2(k);  

    R2_2 = XY1.^2 + XY2.^2 + ep2^2; 
    R_2 = sqrt( R2_2 ); 

    [H1_2, H2_2, S1_2, S2_2] = reg_fncs_withdoublet(ep2, R_2, blob_num); %Mx1

 
    norm1=normal(k,1); 
    norm2=normal(k,2); 
    
    normxy=norm1*XY1+norm2*XY2;

    u1(:)=u1(:)-beta(k)*norm1*(S1_2.*norm1+S2_2.*normxy.*XY1)*g1(k)+...
        -beta(k)*norm2*(S2_2.*normxy.*XY1)*g2(k);

    u2(:)=u2(:)-beta(k)*norm2*(S1_2.*norm2+S2_2.*normxy.*XY2)*g2(k)+...
        -beta(k)*norm1*(S2_2.*normxy.*XY2)*g1(k);    


end




u1 = u1/(mu); 
u2 = u2/(mu); 

%repacking output 
u_beta = [u1 u2]; 



