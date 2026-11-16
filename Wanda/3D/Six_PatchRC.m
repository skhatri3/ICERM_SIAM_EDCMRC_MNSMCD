function [Px,Py,Pz,dArea]=Six_PatchRC(h,A,B)

% this routine creates a discretization on the surface 
% of the unit sphere using the six-patch structure:
%
% Input: 
%      h = the grid space on each side of the cube
%
%When we discretize the cube to find the sphere discretization,
%it is better to have the points on the faces of the cube at the
%center of the grid cells.
%    A = -1;  B = 1;
    
    x = (A+h/2 : h : B-h/2)';
    [a1,a2] = meshgrid(x,x); a3 = ones(size(a1));
    
    am = a1(:)-h/2; bm = a2(:)-h/2;
    ap = a1(:)+h/2; bp = a2(:)+h/2;
    
    %COMPUTE THE EXACT SURFACE AREA OF EACH PROJECTED GRID CELL
    FF1 = 1/2*atan2(2*ap.*bp.*sqrt(1+ap.^2+bp.^2),1+bp.^2-ap.^2.*(-1+bp.^2));
    FF2 = 1/2*atan2(2*am.*bp.*sqrt(1+am.^2+bp.^2),1+bp.^2-am.^2.*(-1+bp.^2));
    FF3 = 1/2*atan2(2*ap.*bm.*sqrt(1+ap.^2+bm.^2),1+bm.^2-ap.^2.*(-1+bm.^2));
    FF4 = 1/2*atan2(2*am.*bm.*sqrt(1+am.^2+bm.^2),1+bm.^2-am.^2.*(-1+bm.^2));
    ddArea = (FF1-FF2)-(FF3-FF4);
    %------------------
    
    xb=[a3(:);-a3(:); a1(:);a1(:); a1(:);a1(:)];
    yb=[a1(:); a1(:);-a3(:);a3(:); a2(:);a2(:)];
    zb=[a2(:); a2(:); a2(:);a2(:);-a3(:);a3(:)];
    dArea = [ddArea;ddArea;ddArea;ddArea;ddArea;ddArea];
    
    bmag = sqrt( xb.^2 + yb.^2 + zb.^2 );
    Px = (xb./bmag);
    Py = (yb./bmag);
    Pz = (zb./bmag);

end