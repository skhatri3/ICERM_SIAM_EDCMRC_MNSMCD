function [errorTable, errorMat] = calculate_errors_fun(U_approx, U_exact, P_approx, P_exact, dx,dy)

% Computes LInf, L1, and L2 errors in velocity and pressure

% Unpack vectors
u_approx = U_approx(:,1);
v_approx = U_approx(:,2);
u_exact = U_exact(:,1);
v_exact = U_exact(:,2);

% Differences
u_error = u_approx - u_exact;
v_error = v_approx - v_exact;
p_error = P_approx - P_exact;

% LInf
u_error_max = max(max(abs(u_error)));
v_error_max = max(max(abs(v_error)));
p_error_max = max(max(abs(p_error)));
% L1
u_error_L1 = dx*dy*sum(sum(abs(u_error)));
v_error_L1 = dx*dy*sum(sum(abs(v_error)));
p_error_L1 = dx*dy*sum(sum(abs(p_error)));
% L2
u_error_L2 = sqrt(dx*dy*sum(sum(u_error.^2)));
v_error_L2 = sqrt(dx*dy*sum(sum(v_error.^2)));
p_error_L2 = sqrt(dx*dy*sum(sum(p_error.^2)));

% organize error calculations
ColNames = {'Max','L1','L2'};
Emax = [u_error_max; v_error_max; p_error_max];
EL1 = [u_error_L1; v_error_L1; p_error_L1];
EL2 = [u_error_L2; v_error_L2; p_error_L2];
errorTable = table(Emax, EL1, EL2,'VariableNames',ColNames);
disp(errorTable)
errorMat = [Emax, EL1, EL2];

end