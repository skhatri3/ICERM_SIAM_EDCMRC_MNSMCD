function [p] = RegStokeslets2D_forcetopressure(y,f,x,ep,blob_num, b, normal)

% Computes pressures at set of points when given a set of points and
% forces at those points using the Method of Regularized Stokeslets 
% Based on Cortez, Fluids, 2021
% we are assuming the source points are an ellipse
a = max(y(:,1));
b = max(y(:,2));

% Force, normals, b(s) to pressure 
% Wanda Strychalski July 31, 2024

%y = (y1,y2) source points of force
%f = (f1,f2) forces at those source points (not force density)
%x = (x1,x2) target points 
%p = pressure evaluated at those target points 
%mu is the viscosity 
%ep is the width of the regularization 
% b(s) is related to the permeability
% normal constains the normal vectors

N = size(y,1); %number of source points 
M = size(x,1); %number of target points

test = x(:,1).^2./a^2 + x(:,2).^2./b^2;
inside = (test<=1);




%unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2);
f1 = f(:,1);
f2 = f(:,2); 
x1 = x(:,1); 
x2 = x(:,2); 
n1 = normal(:,1);
n2 = normal(:,2);
fn = dot(f,normal,2); % check the size

% target points
jumpidx = zeros(size(x));
[ig,jg] = size(y);
xf=y(:,1);   yf=y(:,2); %NL=fix(length(yf)*4/11);
for k1 = 1 : ig
    for k2 = 1: jg
        tmp = sqrt((xg(k1,k2)-xf(1:2*NL)).^2+(yg(k1,k2)-yf(1:2*NL)).^2);
         jk = find( tmp==min(tmp) , 1 );
         jumpidx(k1,k2) = jk;
    end
end

%initializing the velocity 
p = zeros(M,1);

%loop over source points    
for k = 1:N 

    %distance between target and source points 
    XY1 = x1(:) - y1(k); 
    XY2 = x2(:) - y2(k);  
    R2 = XY1.^2 + XY2.^2 + ep^2; 
    R = sqrt( R2 ); 

    %computing the pressure
    %Note: S1 = G'(r)/r
    %p(r) = (f . x)S1
    [~, ~, S1, ~, ~] = reg_fncs_withdoublet(ep,R,blob_num);

    fdotXY = f1(k)*XY1 + f2(k)*XY2; 
    
    p(:) = p(:) + fdotXY.*S1;

end

   pg(k1,k2) = sum( fdotx.*Hs/mu )*hL ...
               + flag*fn(jumpidx(k1,k2))/1.*(yg(k1,k2)>0).*(yg(k1,k2)<H);

end