function [xi, eta, w] = getTriangleQuadrature()
    % 7-point Dunavant Rule (Degree of precision = 5)
    % Defined on the reference triangle: 0 <= xi, eta <= 1, xi + eta <= 1
    
    % Barycentric coordinates (lambda1, lambda2, lambda3) for the 7 points
    % Point 1: Center point
    c1 = 1/3;
    lambda(1,:) = [c1, c1, c1];
    w(1) = 0.225000000000000;
    
    % Points 2,3,4: First symmetric group
    a1 = 0.059715871789770; b1 = 0.470142064105115;
    lambda(2,:) = [a1, b1, b1];
    lambda(3,:) = [b1, a1, b1];
    lambda(4,:) = [b1, b1, a1];
    w(2:4) = 0.132394152788476;
    
    % Points 5,6,7: Second symmetric group
    a2 = 0.797426985353087; b2 = 0.101286507323456;
    lambda(5,:) = [a2, b2, b2];
    lambda(6,:) = [b2, a2, b2];
    lambda(7,:) = [b2, b2, a2];
    w(5:7) = 0.125939180544124;
    
    % Convert Barycentric (lambda1, lambda2) to Cartesian reference (xi, eta)
    xi  = lambda(:, 1);
    eta = lambda(:, 2);
    w   = w(:); % Ensure column vector
end