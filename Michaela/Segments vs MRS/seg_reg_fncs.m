function [H10, H11, H20, H21, H22, H23, S10, S11, S20, S21, S22] = seg_reg_fncs(x,y,ep,blob_num)

%   * x = [x1, x2] are the evaluation/target points
%   * y = [yj,yk] are the segment end points
%   * ep is the epsilon value for the segment
%   * blob_num is the blob choice, 1 = phi (basic), 2 = psi (fancy)
%
% Needs access to file segment_int_terms.m
%
% Outputs the components needed for the segment method solution

ep2 = ep.^2;

switch blob_num
    case 1 % The basic blob
         % phi = (2*d^4)/(pi*(r^2+d^2)^3) 
         [Q0,Q1,T02,T04,T06,~,~,T12,T14,T16,~,~,T22,T24,T26,~,~,T32,~,~] = segment_int_terms(x,y,ep);
         H10 = (ep2.*T02 - Q0)/(4*pi);
         H11 = (ep2.*T12 - Q1)/(4*pi);
         H20 = (T02)/(4*pi);
         H21 = (T12)/(4*pi);
         H22 = (T22)/(4*pi);
         H23 = (T32)/(4*pi);
         S10 = (T02 + ep2.*T04)/(2*pi);
         S11 = (T12 + ep2.*T14)/(2*pi);
         S20 = -(4*ep2.*T06 + 2*T04)/(2*pi);
         S21 = -(4*ep2.*T16 + 2*T14)/(2*pi);
         S22 = -(4*ep2.*T26 + 2*T24)/(2*pi);
    case 2 % The fancy blob
        % psi = 2*ep^4*(r^4 - 10*ep^2*r^2 + 5*ep^4)/(Pi*(r^2 + ep^2)^5)
        [Q0,Q1,T02,T04,T06,T08,T010,T12,T14,T16,T18,T110,T22,T24,T26,~,T210,T32,T34,T36] = segment_int_terms(x,y,ep);
        ep4 = ep2.*ep2; ep6 = ep4.*ep2;
        H10 = ((1/3)*ep2.*T02 - (2/3)*ep4.*T04 + (8/3)*ep6.*T06 - Q0)/(4*pi);
        H11 = ((1/3)*ep2.*T12 - (2/3)*ep4.*T14 + (8/3)*ep6.*T16 - Q1)/(4*pi);
        H20 = (T02 + (4/3)*ep2.*T04 + (8/3)*ep4.*T06)/(4*pi);
        H21 = (T12 + (4/3)*ep2.*T14 + (8/3)*ep4.*T16)/(4*pi);
        H22 = (T22 + (4/3)*ep2.*T24 + (8/3)*ep4.*T26)/(4*pi);
        H23 = (T32 + (4/3)*ep2.*T34 + (8/3)*ep4.*T36)/(4*pi);
        S10 = (T02 + ep2.*T04 + 8*ep6.*T08)/(2*pi);
        S11 = (T12 + ep2.*T14 + 8*ep6.*T18)/(2*pi);
        S20 = -(T04 + 2*ep2.*T06 + 32*ep6.*T010)/pi;
        S21 = -(T14 + 2*ep2.*T16 + 32*ep6.*T110)/pi;
        S22 = -(T24 + 2*ep2.*T26 + 32*ep6.*T210)/pi;
end