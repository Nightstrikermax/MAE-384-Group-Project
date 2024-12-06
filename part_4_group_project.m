%%Problem 4
% Parameters
Total_pop = 1000; % Total population
Sim_days = 30;   % Duration of the simulation in days
T_step = 0.1;       % Time step in days

% Initial conditions
S0 = 990; %initial suceptible
I0 = 10;  %initial infected
R0 = 0;   %initial recovered

%Trans rate and recovery rate
trans_rate = @ (x) (0.3*(1+5*sin(2*pi*x))); %transmission raten1
rec_rate = 0.1;   %recovery rate

% Initialize arrays for S, I, R values
susceptible = zeros(1, Sim_days/T_step + 1);
infected = zeros(1, Sim_days/T_step + 1);
recovered = zeros(1, Sim_days/T_step + 1);

% Set initial values
susceptible(1) = S0;
infected(1) = I0;
recovered(1) = R0;

n = 1;

% Using 4th-order Runge-Kutta method
for t = 0:0.1:Sim_days-0.1
    % Calculate k1 values
    k1_s = -trans_rate(t) * susceptible(n) * infected(n) / Total_pop;
    k1_i = (trans_rate(t) * susceptible(n) * infected(n) / Total_pop) - rec_rate * infected(n);
    k1_r = rec_rate * infected(n);

    % Calculate k2 values
    k2_s = -trans_rate(t) * (susceptible(n) + 0.5 * T_step * k1_s) * (infected(n) + 0.5 * T_step * k1_i) / Total_pop;
    k2_i = (trans_rate(t) * (susceptible(n) + 0.5 * T_step * k1_s) * (infected(n) + 0.5 * T_step * k1_i) / Total_pop) ...
        - rec_rate * (infected(n) + 0.5 * T_step * k1_i);
    k2_r = rec_rate * (infected(n) + 0.5 * T_step * k1_i);

    % Calculate k3 values
    k3_s = -trans_rate(t) * (susceptible(n) + 0.5 * T_step * k2_s) * (infected(n) + 0.5 * T_step * k2_i) / Total_pop;
    k3_i = (trans_rate(t) * (susceptible(n) + 0.5 * T_step * k2_s) * (infected(n) + 0.5 * T_step * k2_i) / Total_pop) ...
        - rec_rate * (infected(n) + 0.5 * T_step * k2_i);
    k3_r = rec_rate * (infected(n) + 0.5 * T_step * k2_i);

    % Calculate k4 values
    k4_s = -trans_rate(t) * (susceptible(n) + T_step * k3_s) * (infected(n) + T_step * k3_i) / Total_pop;
    k4_i = (trans_rate(t) * (susceptible(n) + T_step * k3_s) * (infected(n) + T_step * k3_i) / Total_pop) ...
        - rec_rate * (infected(n) + T_step * k3_i);
    k4_r = rec_rate * (infected(n) + T_step * k3_i);

    % Update values using weighted average of k1, k2, k3, k4
    susceptible(n + 1) = susceptible(n) + (T_step / 6) * (k1_s + 2 * k2_s + 2 * k3_s + k4_s);
    infected(n + 1) = infected(n) + (T_step / 6) * (k1_i + 2 * k2_i + 2 * k3_i + k4_i);
    recovered(n + 1) = recovered(n) + (T_step / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r);

    n=n+1;
end

% Generate plot
figure;
plot(0:0.1:Sim_days, susceptible, 'b-', 'LineWidth', 1.5); hold on;
plot(0:0.1:Sim_days, infected, 'r-', 'LineWidth', 1.5);
plot(0:0.1:Sim_days, recovered, 'g-', 'LineWidth', 1.5);
hold off;
xlabel('Days');
ylabel('Population');
legend('Susceptible', 'Infected', 'Recovered');
title(['SIR Model (\beta = ', 'Periodic variation', ', \gamma = ', num2str(rec_rate), ')']);


% Fourier transform
susceptiblefft = fft(susceptible);
infectedfft = fft(infected);
recoveredfft = fft(recovered);

% Define frequency
T = Sim_days;
N = 300;
f = 1/T.*(0:N/2);

% Plot
figure;
plot(f,abs(infectedfft(1:N/2+1)))
hold on
title('Spectrum with ω = 2pi')

%% Replace ω with ω = 2π × 100/365

%%Problem 4
% Parameters
Total_pop = 1000; % Total population
Sim_days = 30;   % Duration of the simulation in days
T_step = 0.1;       % Time step in days

% Initial conditions
S0 = 990; %initial suceptible
I0 = 10;  %initial infected
R0 = 0;   %initial recovered

%Trans rate and recovery rate
trans_rate = @ (x) (0.3*(1+5*sin(2*100/365*pi*x))); %transmission raten1
rec_rate = 0.1;   %recovery rate

% Initialize arrays for S, I, R values
susceptible1 = zeros(1, Sim_days/T_step + 1);
infected1 = zeros(1, Sim_days/T_step + 1);
recovered1 = zeros(1, Sim_days/T_step + 1);

% Set initial values
susceptible1(1) = S0;
infected1(1) = I0;
recovered1(1) = R0;

n = 1;

% Using 4th-order Runge-Kutta method
for t = 0:0.1:Sim_days-0.1
    % Calculate k1 values
    k1_s = -trans_rate(t) * susceptible1(n) * infected1(n) / Total_pop;
    k1_i = (trans_rate(t) * susceptible1(n) * infected1(n) / Total_pop) - rec_rate * infected1(n);
    k1_r = rec_rate * infected1(n);

    % Calculate k2 values
    k2_s = -trans_rate(t) * (susceptible1(n) + 0.5 * T_step * k1_s) * (infected1(n) + 0.5 * T_step * k1_i) / Total_pop;
    k2_i = (trans_rate(t) * (susceptible1(n) + 0.5 * T_step * k1_s) * (infected1(n) + 0.5 * T_step * k1_i) / Total_pop) ...
        - rec_rate * (infected1(n) + 0.5 * T_step * k1_i);
    k2_r = rec_rate * (infected1(n) + 0.5 * T_step * k1_i);

    % Calculate k3 values
    k3_s = -trans_rate(t) * (susceptible1(n) + 0.5 * T_step * k2_s) * (infected1(n) + 0.5 * T_step * k2_i) / Total_pop;
    k3_i = (trans_rate(t) * (susceptible1(n) + 0.5 * T_step * k2_s) * (infected1(n) + 0.5 * T_step * k2_i) / Total_pop) ...
        - rec_rate * (infected1(n) + 0.5 * T_step * k2_i);
    k3_r = rec_rate * (infected1(n) + 0.5 * T_step * k2_i);

    % Calculate k4 values
    k4_s = -trans_rate(t) * (susceptible1(n) + T_step * k3_s) * (infected1(n) + T_step * k3_i) / Total_pop;
    k4_i = (trans_rate(t) * (susceptible1(n) + T_step * k3_s) * (infected1(n) + T_step * k3_i) / Total_pop) ...
        - rec_rate * (infected1(n) + T_step * k3_i);
    k4_r = rec_rate * (infected1(n) + T_step * k3_i);

    % Update values using weighted average of k1, k2, k3, k4
    susceptible1(n + 1) = susceptible1(n) + (T_step / 6) * (k1_s + 2 * k2_s + 2 * k3_s + k4_s);
    infected1(n + 1) = infected1(n) + (T_step / 6) * (k1_i + 2 * k2_i + 2 * k3_i + k4_i);
    recovered1(n + 1) = recovered1(n) + (T_step / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r);

    n=n+1;
end

% Generate plot
figure;
plot(0:0.1:Sim_days, susceptible1, 'b-', 'LineWidth', 1.5); hold on;
plot(0:0.1:Sim_days, infected1, 'r-', 'LineWidth', 1.5);
plot(0:0.1:Sim_days, recovered1, 'g-', 'LineWidth', 1.5);
hold off;
xlabel('Days');
ylabel('Population');
legend('Susceptible', 'Infected', 'Recovered');
title(['SIR Model (\beta = ', 'Periodic variation', ', \gamma = ', num2str(rec_rate), ')']);


% Fourier transform
susceptiblefft1 = fft(susceptible1);
infectedfft1 = fft(infected1);
recoveredfft1 = fft(recovered1);

% Define frequency
T = Sim_days;
N = 300;
f = 1/T.*(0:N/2);

% Plot
figure;
plot(f,abs(infectedfft(1:N/2+1)))
hold on
title('Spectrum with ω = 2pi*100/365')


