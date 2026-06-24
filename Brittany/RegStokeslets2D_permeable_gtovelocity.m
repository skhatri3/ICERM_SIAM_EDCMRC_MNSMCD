function [u_beta] = RegStokeslets2D_permeable_gtovelocity(y,g,...
    x,ep,mu,blob_num, beta, normal)

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
%ep is the width of the regularization 
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


%Build Matrix
Beta=zeros(M,N);
Norm1=zeros(M,N);
Norm2=zeros(M,N);
%loop over source points    
for k = 1:N 

    %distance between target and source points 
    XY1(:,k) = x1(:) - y1(k); 
    XY2(:,k) = x2(:) - y2(k);  

    %Stuff needed for doublet:
    %Beta matrix: column k has beta_k
    Beta(:, k)=beta(k);

    %Normal matrices needed: 
    Norm1(:,k)=normal(k,1); %kth column is all nk1
    Norm2(:,k)=normal(k,2); %kth column is all nk2
    
   
end



R2 = XY1.^2 + XY2.^2 + ep^2; 
R = sqrt( R2 ); 

[H1, H2, S1, S2] = reg_fncs_withdoublet(ep, R, blob_num); %MxN


%Doublet part
% S1(S1(:)==-Inf|S1(:)==Inf)=0;
% S2(S2(:)==-Inf|S2(:)==Inf)=0;
NormXY=Norm1.*XY1+Norm2.*XY2;
D11=-Beta.*Norm1.*(S1.*Norm1+S2.*NormXY.*XY1);
% D12=-Beta.*Norm2.*(S1.*Norm1+S2.*NormXY.*XY1);
% D21=-Beta.*Norm1.*(S1.*Norm2+S2.*NormXY.*XY2);
D12=-Beta.*Norm2.*(S2.*NormXY.*XY1);
D21=-Beta.*Norm1.*(S2.*NormXY.*XY2);

D22=-Beta.*Norm2.*(S1.*Norm2+S2.*NormXY.*XY2);
D=[D11 D12; D21 D22]/(mu);
% %zero out rows for targets not on permeable membrane
% D(I2,:)=zeros(length(I2),length(D(1,:)));
% D(I2+M,:)=zeros(length(I2), length(D(1,:)));

%Obtain velocities
gg=[g1; g2];
uu=D*gg;
%repack
u1=uu(1:M);
u2=uu(M+1:end);
u_beta=[u1,u2];


