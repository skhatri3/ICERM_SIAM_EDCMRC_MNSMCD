function [Q0,Q1,T02,T04,T06,T08,T010,T12,T14,T16,T18,T110,T22,T24,T26,T28,T210,T32,T34,T36] = segment_int_terms(x,y,ep)
%
%   * x = [x1, x2] are the evaluation/target points
%   * y = [yj,yk] are the segment end points
%   * ep is the epsilon value for the segment
% 
% Inputs x and y can be ( )x2 matrices, so that you can output T's for
% multiple segments and evaluation points. ep should be the same length as
% y.
%
% Developed by Michaela Kubacki and Ricardo Cortez March 2026

yj = y(1,:);
yk = y(2,:);
x1 = x(:,1);
x2 = x(:,2);

xj1 = x1 - yj(1);
xj2 = x2 - yj(2);
xk1 = x1 - yk(1);
xk2 = x2 - yk(2);

ep2 = ep.^2;

% R^2 = r^2 + ep^2, below R is evaluated at alpha = 0 and alpha = 1
R0 = sqrt(xj1.^2 + xj2.^2 + ep2);
R1 = sqrt(xk1.^2 + xk2.^2 + ep2);

v = yj-yk; % Tangent vector
L = norm(v); % Segment length
L2 = L^2;

xjdotv = xj1*v(1) + xj2*v(2);
xkdotv = xk1*v(1) + xk2*v(2);

Log0 = log(R0);
Log1 = log(R1);

c = R0.^2 - (xjdotv).^2/L2;
csqrt = sqrt(c);

% Change of variables in formulas: xi = alpha*L + (xjdotv)/L

% Evaluating xi at alpha = 0 and alpha = 1
xi1 = L + (xjdotv)/L;
xi0 = xjdotv/L;

% Evaluating arctan solution component when alpha = 0 and alpha = 1
Atan1 = atan2(xi1,csqrt);
Atan0 = atan2(xi0,csqrt);

% Powers of R evaluated at alpha = 0 and alpha = 1
R02 = R0.^2; R04 = R02.*R02; R06 = R02.*R04; R08 = R02.*R06;
R12 = R1.^2; R14 = R12.*R12; R16 = R12.*R14; R18 = R12.*R16;

% Q0 calculated directly (1) in overleaf
Q0=1/L2*( xkdotv.*(Log1-1) - xjdotv.*(Log0-1) + L*csqrt.*(Atan1-Atan0) );

% T02 calculated directly (1) in overleaf
T02 = (Atan1-Atan0)./(L*csqrt);

% Q1 using recursion with Q0 (4) in overleaf
Q1 = (2*R12.*Log1 - R12)/(4*L2) - (2*R02.*Log0 - R02)/(4*L2) - 1/(L2)*xjdotv.*Q0;

% Rest of T0's calculated using recursion T0(2n) recursion (2) in overleaf
T04 = (0.5*T02 + 0.5*xi1./(R12)/L - 0.5*xi0./(R02)/L)./c;
T06 = (0.75*T04 + 0.25*xi1./(R14)/L - 0.25*xi0./(R04)/L )./c;
T08 = ((5/6)*T06 + (1/6)*xi1./R16/L - (1/6)*xi0./R06/L)./c;
T010 = ((7/8)*T08 + (1/8)*xi1./R18/L - (1/8)*xi0./R08/L)./c;

% Rest of T's Calculated using the recursion (5) in overleaf
T12 = (Log1-Log0)/L2 - xjdotv.*T02/L2;
T22 = Log1/L2 - Q0/L2 - xjdotv.*T12/L2;
T32 = Log1/L2 - 2*Q1/L2 - xjdotv.*T22/L2;

T14 = -0.5/L2 * (1./R12 - 1./R02) - xjdotv.*T04/L2;
T24 = -0.5/L2 * (1./R12) + 0.5/L2 * T02 - xjdotv.*T14/L2;
T34 = -0.5/L2 * (1./R12) + T12/L2 - xjdotv.*T24/L2;

T16 = -0.25/L2 * (1./R14 - 1./R04) - xjdotv.*T06/L2;
T26 = -0.25/L2 * (1./R14) + 0.25/L2 * T04 - xjdotv.*T16/L2;
T36 = -0.25/L2 * (1./R14) + 0.5/L2 * T14 - xjdotv.*T26/L2;

T18 = -1/6/L2 * (1./R16 - 1./R06) - xjdotv.*T08/L2;
T28 = -1/6/L2 * (1./R16) + (1/6)/L2 * T06 - xjdotv.*T18/L2;

T110 = -1/8/L2 * (1./R18 - 1./R08) - xjdotv.*T010/L2;
T210 = -1/8/L2 * (1./R18) + (1/8)/L2 *T08 - xjdotv.*T110/L2;



