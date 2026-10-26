function [errorTable, errorMat] = calculate_errors_fun(U_approx, U_exact, P_approx, P_exact, dx, dy)

% Computes LInf, L1, and L2 errors in velocity and pressure

haveU = ~(isempty(U_approx) || isempty(U_exact));
haveP = ~(isempty(P_approx) || isempty(P_exact));

% initialize vectors for storing error calculations
Emax = [];
EL1 = [];
EL2 = [];
RowNames = {};

if ~haveU && ~haveP
    error('Provide at least one of the U or P pairs.');
end

if haveU
    % Unpack vectors
    u_approx = U_approx(:,1);
    v_approx = U_approx(:,2);
    u_exact = U_exact(:,1);
    v_exact = U_exact(:,2);

    % Differences
    u_error = u_approx - u_exact;
    v_error = v_approx - v_exact;

    % LInf
    u_error_max = max(max(abs(u_error)));
    v_error_max = max(max(abs(v_error)));
    
    % L1
    u_error_L1 = sum(sum(abs(u_error)))/(Nx1*Nx2); %dx*dy*sum(sum(abs(u_error)));
    v_error_L1 = sum(sum(abs(v_error)))/(Nx1*Nx2); %dx*dy*sum(sum(abs(v_error)));
    
    % L2
    u_error_L2 = sqrt(sum(sum(u_error.^2))/(Nx1*Nx2)); %sqrt(dx*dy*sum(sum(u_error.^2)));
    v_error_L2 = sqrt(sum(sum(v_error.^2))/(Nx1*Nx2)); %sqrt(dx*dy*sum(sum(v_error.^2)));

    % organize error calculations
    Emax = [Emax; u_error_max; v_error_max];
    EL1 = [EL1; u_error_L1; v_error_L1];
    EL2 = [EL2; u_error_L2; v_error_L2];
    RowNames = [RowNames, {'u'}, {'v'}];
end

if haveP
    % Difference
    p_error = P_approx - P_exact;

    % LInf
    p_error_max = max(max(abs(p_error)));

    % L1
    p_error_L1 = dx*dy*sum(sum(abs(p_error)));
    
    % L2
    p_error_L2 = sqrt(dx*dy*sum(sum(p_error.^2)));

    % organize error calculations
    Emax = [Emax; p_error_max];
    EL1 = [EL1; p_error_L1];
    EL2 = [EL2; p_error_L2];

    RowNames = [RowNames, {'p'}];
end

% Construct error table and matrix
ColNames = {'Max','L1','L2'};
errorTable = table(Emax, EL1, EL2, 'VariableNames', ColNames, 'RowNames', RowNames);
disp(errorTable)
errorMat = [Emax, EL1, EL2];

end