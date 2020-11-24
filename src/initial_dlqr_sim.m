%%% Inverted Pendulum DLQR Simulation

% Parameters (SI)
m = 1; 
M = 20;
L = 0.5;
g = 9.81;

% State-Space Matrices
A = [0 1 0 0;
     0 0 m*g/M 0;
     0 0 0 1;
     0 0 (M+m)*g/(M*L) 0];
 
B = [0 1/M 0 1/(M*L)];

