# MAE 384 Group Project
 This github repository is for Matlab code which models the spraed of disease using the SIR model. This stands for Susceptible-Infected-Recovered which is used for many modern day epidemics. The goal of this code is to autonomously project the amount of people affected and how 
Part 1:
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


Part 4:
Part 4 was done in mostly the same way as part 1 but instead of the transmition rate being constant a function handle was used to find and used the transmition value for each time the loop was ran. Part 4 also included MATLABS fft function to plot the spectrum. The fft (Fast Fourier Transofrm) function converts discrete signals from the time domain to the frequency domain and provides information about the frequency content, phase, and other properties of the signal.
