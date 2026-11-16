function [q] = RegStokeslets2D_permeable_forcetopressure(y,x,ep,blob_num,...
    b,normal)

% Computes permeable pressure q
% Modified by Wanda Strychalski July 31, 2024

% y = (y1,y2) source points of force
% x = (x1,x2) target points
% q = pressure evaluated at those target points
% ep is the width of the regularization
% b is b(s) in the paper
% n is Nx2, giving unit normals for source points

N = size(y,1); %number of source points
M = size(x,1); %number of target points

% for k1 = 1 : N
%     for k2 = 1: M
%         tmp = sqrt((xg(k1,k2)-xf(1:2*NL)).^2+(yg(k1,k2)-yf(1:2*NL)).^2);
%          jk = find( tmp==min(tmp) , 1 );
%          jumpidx(k1,k2) = jk;
%     end
% end

%unpacking the inputs
y1 = y(:,1);
y2 = y(:,2);
% f1 = f(:,1);
% f2 = f(:,2);
x1 = x(:,1);
x2 = x(:,2);
n1 = normal(:,1);
n2 = normal(:,2);

%initializing the velocity
q = zeros(M,1);

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
    [~, ~, S1, ~,Q] = reg_fncs_withdoublet(ep,R,blob_num);

    % fdotXY = f1(k)*XY1 + f2(k)*XY2;
    ndotXY = n1(k)*XY1+n2(k)*XY2;

    q(:) = q(:) - b.*ndotXY.*Q;

end

end