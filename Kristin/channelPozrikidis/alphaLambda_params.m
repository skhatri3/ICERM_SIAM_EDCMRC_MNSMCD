% Check relationship between alpha and lambda for Da --> 0

% Darcy number
Da_vec = 0.001:0.001:0.1;

alpha_vec = zeros(size(Da_vec));
lambda_vec = zeros(size(Da_vec));
gamma_vec = zeros(size(Da_vec));

% Channel geometry
L = 2; % length of permeable portion of the channel
H = 1; % radius of the channel
c = 2*H; % extension length of the channel
G = 2; % Poiseuille flow strength
mu = 1; % viscosity

for j = 1:length(Da_vec)

    Da = Da_vec(j);
    % Compute alpha and gamma from CP methodology
    alpha = sqrt(3*mu*Da/(H^3)); % needed for CP solution
    gamma = -G/(alpha*cosh(alpha*L)); % needed for CP solution
    alpha_vec(j) = alpha;
    gamma_vec(j) = gamma;

    % Use Newton's method to compute lambda from FB methodology
    l_old = 0.5;
    for nw = 1:30
        f = (Da-H/2)*l_old*cot(l_old*H)*cos(l_old*H) + cos(l_old*H)/2 - H*l_old*sin(l_old*H)/2;
        fp = cos(l_old*H)*((Da-H/2)*(cot(l_old*H)-H*l_old*csc(l_old*H)^2)-H^2*l_old/2) - H*sin(l_old*H)*(l_old*cot(l_old*H)+1);
        lambda = l_old - f/fp;
        l_old = lambda;
    end
    lambda_vec(j) = lambda;

end

%% Plot
figure
plot(Da_vec, lambda_vec, '-', 'LineWidth', 1.3)
hold on
plot(Da_vec, alpha_vec, '-', 'LineWidth', 1.3)
xlabel('Da')
legend('$\lambda$', '$\alpha$', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', 20)
ax = gca; ax.FontSize = 16;

paramDiffs = abs(alpha_vec - lambda_vec);
figure
plot(Da_vec, paramDiffs, '-', 'LineWidth', 1.3)
xlabel('Da')
ylabel('$|\alpha-\lambda|$', 'Interpreter', 'latex')
ax = gca; ax.FontSize = 16;