%%MAE384 Group Project PART II: Interpolation


beta = 0.3;  % rate of infection 
gamma = 0.1; % rate of recovery 
N = 1;       % size of population 
T = 100;     % total simulation 
h1 = 1;      % finer
h2 = 2;      % coarser 


S0 = 0.99; % susceptible 
I0 = 0.01; % infected
R0 = 0;    % recovered 


t_f = 0:h1:T; % fine steps
t_c = 0:h2:T; % coarse steps

%  time step (finer)
S_f = zeros(size(t_f));
I_f = zeros(size(t_f));
R_f = zeros(size(t_f));

% initial conditions 
S_f(1) = S0;
I_f(1) = I0;
R_f(1) = R0;

% time step 
for k = 1:length(t_fine)-1
    dS = -(beta/N) * S_f(k) * I_f(k);
    dI = (beta/N) * S_f(k) * I_f(k) - gamma * I_f(k);
    dR = gamma * I_f(k);
    
    S_f(k+1) = S_f(k) + h1 * dS;
    I_f(k+1) = I_f(k) + h1 * dI;
    R_f(k+1) = R_f(k) + h1 * dR;
end

% coarser step
S_c = zeros(size(t_c));
I_c= zeros(size(t_c));
R_c = zeros(size(t_c));

%intial conditons 
S_c(1) = S0;
I_c(1) = I0;
R_c(1) = R0;

% time step
for k = 1:length(t_c)-1
    dS = -(beta/N) * S_c(k) * I_c(k);
    dI = (beta/N) * S_c(k) * I_c(k) - gamma * I_c(k);
    dR = gamma * I_c(k);
    
    S_c(k+1) = S_c(k) + h2 * dS;
    I_c(k+1) = I_c(k) + h2 * dI;
    R_c(k+1) = R_c(k) + h2 * dR;
end

% interpolation of odd days
t_odd = 1:2:T-1;

% coaser linear interpolation
S_l = interp1(t_c, S_c, t_odd, 'linear');
I_l = interp1(t_c, I_c, t_odd, 'linear');
R_l = interp1(t_c, R_c, t_odd, 'linear');

% lagrange 
S_q = interp1(t_c, S_coarser, t_odd, 'spline');
I_q = interp1(t_c, I_c, t_odd, 'spline');
R_q = interp1(t_c, R_c, t_odd, 'spline');

% finer odd
S_f_odd = interp1(t_f, S_f, t_odd);
I_f_odd = interp1(t_f, I_f, t_odd);
R_f_odd = interp1(t_f, R_f, t_odd);

% linear interpolation
Nint = length(t_odd);
EL2_S_l = sqrt(sum((S_l - S_f_odd).^2) / Nint);
EL2_I_l = sqrt(sum((I_l - I_f_odd).^2) / Nint);
EL2_R_l = sqrt(sum((R_l - R_f_odd).^2) / Nint);

% quad. interpolation
el2_S_q = sqrt(sum((S_q - S_f_odd).^2) / Nint);
el2_I_q = sqrt(sum((I_q - I_f_odd).^2) / Nint);
el2_R_q = sqrt(sum((R_q - R_f_odd).^2) / Nint);

% error table
ErrorTable = table(["Linear"; "Quadratic"], ...
                   [EL2_S_l; el2_S_q], ...
                   [EL2_I_l; el2_I_q], ...
                   [EL2_R_l; el2_R_q], ...
                   'VariableNames', {'Interpolation', 'S_Error', 'I_Error', 'R_Error'});

% error table
disp(ErrorTable);

% plot 
figure;
plot(t_odd, S_f_odd, 'k-o', 'DisplayName', 'Fine Solution (Odd Days)');
hold on;
plot(t_odd, S_l, 'b-*', 'DisplayName', 'Linear Interpolation');
plot(t_odd, S_q, 'r-d', 'DisplayName', 'Quadratic Interpolation');
xlabel('Time (days)');
ylabel('S(t)');
legend;
title('Comparison of Interpolations for S(t)');
grid on;
