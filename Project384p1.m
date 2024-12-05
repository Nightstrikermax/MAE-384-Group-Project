%%Problem 1
% Parameters
Total_pop = 1000; % Total population
Sim_days = 100;   % Duration of the simulation in days
T_step = 1;       % Time step in days

% Initial conditions
S0 = 990; %initial suceptible
I0 = 10;  %initial infected
R0 = 0;   %initial recovered

% Parameters for diseases
dis_params = [0.3, 0.1;  % Seasonal Influenza
              1.0, 0.1;  % COVID-19
             2.0, 0.2]; % Measles

% Loop through each disease scenario
for cdc = 1:size(dis_params, 1)
 
    trans_rate = dis_params(cdc, 1); %transmission rate
    rec_rate = dis_params(cdc, 2);   %recovery rate
    
    % Initialize arrays for S, I, R values
    susceptible = zeros(1, Sim_days + 1);
    infected = zeros(1, Sim_days + 1);
    recovered = zeros(1, Sim_days + 1);
    
    % Set initial values
    susceptible(1) = S0;
    infected(1) = I0;
    recovered(1) = R0;
    
    % Using 4th-order Runge-Kutta method
    for t = 1:Sim_days
        % Calculate k1 values
        k1_s = -trans_rate * susceptible(t) * infected(t) / Total_pop;
        k1_i = (trans_rate * susceptible(t) * infected(t) / Total_pop) - rec_rate * infected(t);
        k1_r = rec_rate * infected(t);
        
        % Calculate k2 values
        k2_s = -trans_rate * (susceptible(t) + 0.5 * T_step * k1_s) * (infected(t) + 0.5 * T_step * k1_i) / Total_pop;
        k2_i = (trans_rate * (susceptible(t) + 0.5 * T_step * k1_s) * (infected(t) + 0.5 * T_step * k1_i) / Total_pop) ...
               - rec_rate * (infected(t) + 0.5 * T_step * k1_i);
        k2_r = rec_rate * (infected(t) + 0.5 * T_step * k1_i);
        
        % Calculate k3 values
        k3_s = -trans_rate * (susceptible(t) + 0.5 * T_step * k2_s) * (infected(t) + 0.5 * T_step * k2_i) / Total_pop;
        k3_i = (trans_rate * (susceptible(t) + 0.5 * T_step * k2_s) * (infected(t) + 0.5 * T_step * k2_i) / Total_pop) ...
               - rec_rate * (infected(t) + 0.5 * T_step * k2_i);
        k3_r = rec_rate * (infected(t) + 0.5 * T_step * k2_i);
        
        % Calculate k4 values
        k4_s = -trans_rate * (susceptible(t) + T_step * k3_s) * (infected(t) + T_step * k3_i) / Total_pop;
        k4_i = (trans_rate * (susceptible(t) + T_step * k3_s) * (infected(t) + T_step * k3_i) / Total_pop) ...
               - rec_rate * (infected(t) + T_step * k3_i);
        k4_r = rec_rate * (infected(t) + T_step * k3_i);
        
        % Update values using weighted average of k1, k2, k3, k4
        susceptible(t + 1) = susceptible(t) + (T_step / 6) * (k1_s + 2 * k2_s + 2 * k3_s + k4_s);
        infected(t + 1) = infected(t) + (T_step / 6) * (k1_i + 2 * k2_i + 2 * k3_i + k4_i);
        recovered(t + 1) = recovered(t) + (T_step / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r);
    end
    
    % Generate plots
    figure;
    plot(0:Sim_days, susceptible, 'b-', 'LineWidth', 1.5); hold on;
    plot(0:Sim_days, infected, 'r-', 'LineWidth', 1.5);
    plot(0:Sim_days, recovered, 'g-', 'LineWidth', 1.5);
    hold off;
    xlabel('Days');
    ylabel('Population');
    legend('Susceptible', 'Infected', 'Recovered');
    title(['SIR Model (\beta = ', num2str(trans_rate), ', \gamma = ', num2str(rec_rate), ')']);
end
