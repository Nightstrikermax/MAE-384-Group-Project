# MAE 384 Group Project
 This github repository is for Matlab code which models the spraed of disease using the SIR model. This stands for Susceptible-Infected-Recovered which is used for many modern day epidemics. The goal of this code is to  
 automously project the amount of people affected using transmission rates, recovery rates, total population, number of infected people,and number of people who are susceptible over time


# Part 1
Mathematical model that is described by the equations below
![image](https://github.com/user-attachments/assets/f29d7841-237d-4007-b91a-c03a594f9939)

Code Workflow for part 1:
Initialization:
Defines the total population (𝑁=1000).
Sets initial conditions (𝑆0=990,𝐼0=10,𝑅0=0).
Defines the parameters (𝛽,𝛾) for the three diseases.

Numerical Solver:
Implements the 4th-order Runge-Kutta method to solve the SIR equations iteratively over 100 days.

Visualization:
Generates a plot for each disease showing the time evolution of 𝑆(𝑡),𝐼(𝑡), and𝑅(𝑡).


# Part 2
This section of the code used quadratic interpolation to provide smaller errors compared to linear interpolation. Quadratic interpolation uses higher orders of polynomials to give a closer approximation. Linear interpolation means the function will be linear between two coarse time steps and it performs better if the function is almost linear as well

# Part 3
This part of the code was used to approximate the amount of initially infected individuals using the SIR model with a constant amount of susceptible individuals and the amount of infected at a certain time. This part of the code was also used to help find the transmission rate of a disease based off of the same numbers used to find the individuals who were initially infected 



# Part 4
Part 4 was done in mostly the same way as part 1 but instead of the transmition rate being constant a function handle was used to find and used the transmition value for each time the loop was ran. Part 4 also included MATLABS fft function to plot the spectrum. The fft (Fast Fourier Transofrm) function converts discrete signals from the time domain to the frequency domain and provides information about the frequency content, phase, and other properties of the signal.
