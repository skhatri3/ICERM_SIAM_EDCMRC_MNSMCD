function [u] = RegStokeslets3D_gforceto_velocity_permeable(y,...
    x,g,ep,mu,blob, beta, normal)


% Computes intermediate "g" forces for permeable membrane 
% using the Method of Regularized Stokeslets 
% Based on Cortez, SIAM J. Sci Comput. 2001 and Cortez, Fluids 2021 

% Developed by Shilpa Khatri, Ricardo Cortez, Brittany Leathers, and
% Michaela Kubacki, Wanda Strychalski
% October 2024 

%y = (y1,y2, y3) source points Nx3
%f = (f1,f2, f3) forces at those source points (not force density) 
%x = (x1,x2, x3) target points 
%u = (u1,u2, u3) velocity at those target points 
%mu is the viscosity 
%ep is the width of the regularization 
%I: Indices of x (target) rows that correspond to Gamma_beta, the permeable region
%Beta gives the permeability coefficient for all source points., Nx1
    %beta(i) should be 0 if it is not a permeable section
%n is Nx3, giving unit normals for source points
% WS This only contains the doublet part!


N = size(y,1); %number of source points 
M = size(x,1); %number of target points 

u = zeros(M,3);

%unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2);
y3 = y(:,3);
g1 = g(:,1);
g2 = g(:,2);
g3 = g(:,3);
x1 = x(:,1);
x2 = x(:,2);
x3=  x(:,3);

%building matrix 

Beta=zeros(M,N);
Norm1=zeros(M,N);
Norm2=zeros(M,N);
%loop over source points    
for k = 1:N 

    %distance between target and source points %MxN
    XY1(:,k) = x1(:) - y1(k); 
    XY2(:,k) = x2(:) - y2(k);  
    XY3(:,k) = x3(:) - y3(k);  
    %Stuff needed for doublet:
    %Beta matrix: column k has beta_k
    Beta(:, k)=beta(k);

    %Normal matrices needed: 
    Norm1(:,k)=normal(k,1); %kth column is all nk1
    Norm2(:,k)=normal(k,2); %kth column is all nk2
    Norm3(:,k)=normal(k,3); %kth column is all nk3
    
   
end

R2 = XY1.^2 + XY2.^2 + XY3.^2+ ep^2; 
R = sqrt( R2 ); 

[H1, H2, S1, S2] = reg_fncs_3D(ep, R, blob); %MxN


%Stokeslet part
M11 = H1 + H2.*XY1.*XY1; 
M22 = H1 + H2.*XY2.*XY2; 
M33 = H1 + H2.*XY3.*XY3; 
M12 = H2.*XY1.*XY2; 
M13 = H2.*XY1.*XY3; 
M23 = H2.*XY2.*XY3; 

Mat = [M11 M12 M13; M12 M22 M23; M13 M23 M33]/mu;


%Doublet part
NormXY=Norm1.*XY1+Norm2.*XY2+Norm3.*XY3;
D11= -Beta.*Norm1.*(S1.*Norm1+S2.*NormXY.*XY1 );
D12= -Beta.*Norm2.*(S1.*Norm1+S2.*NormXY.*XY1 );
D13= -Beta.*Norm3.*(S1.*Norm1+S2.*NormXY.*XY1 );
D21= -Beta.*Norm1.*(S1.*Norm2+S2.*NormXY.*XY2 );
D22= -Beta.*Norm2.*(S1.*Norm2+S2.*NormXY.*XY2 );
D23= -Beta.*Norm3.*(S1.*Norm2+S2.*NormXY.*XY2 );
D31= -Beta.*Norm1.*(S1.*Norm3+S2.*NormXY.*XY3 );
D32= -Beta.*Norm2.*(S1.*Norm3+S2.*NormXY.*XY3 );
D33= -Beta.*Norm3.*(S1.*Norm3+S2.*NormXY.*XY3 );

D = [D11 D12 D13; D21 D22 D23; D31 D32 D33]/mu;
%zero out rows that are for target points on permeable membrane
% D(I,:)=zeros(length(I),length(D(1,:)));
% D(I+M,:)=zeros(length(I), length(D(1,:)));
% D(I+2*M,:)=zeros(length(I), length(D(1,:)));

%solving for force when given a velocity 
%uu = [u1;u2;u3];
%gg = (Mat+D)\uu; 
gg = [g1;g2;g3];
%uu = (Mat+D)*gg;
uu = D*gg;% Doublet part only
% g1 = gg(1:N); 
% g2 = gg(N+1:2*N);
% g3 = gg(2*N+1:end);  
u1 = uu(1:N);
u2 = uu(N+1:2*N);
u3 = uu(2*N+1:end);


%repacking output 
g = [g1,g2,g3]; 
%repacking output 
u = [u1,u2,u3]; 
