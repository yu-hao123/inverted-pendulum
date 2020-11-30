%%% Inverted Pendulum DLQR Simulation
clear variables;
close all;

% Parameters (SI)
m = 1; % kg
M = 9; % kg
L = 0.5; % m
g = 9.81; % m/s^2

% State-Space Matrices (Continuous)
Ac = [0 1 0 0;
     0 0 m*g/M 0;
     0 0 0 1;
     0 0 (M+m)*g/(M*L) 0];
 
Bc = [0; 1/M; 0; 1/(M*L)];

Ts = 0.05;
n_samples = 150;

[A, B] = c2dm(Ac, Bc, [], [], Ts, 'zoh');

ubar = 0;
xbar = [1; 0; 0; 0];
x0 = [0; 0; 0; 0];

% Gain Q and R matrices
R = 1 / 10^2;
q1 = 1;
q2 = 1;
q3 = 1 / (pi/10)^2;
q4 = 1 / (pi/10)^2;
Q = diag([q1, q2, q3, q4]);

[K,S,E] = dlqr(A,B,Q,R);

x{1} = x0;
for k = 1:n_samples
    u(k) = ubar - K * (x{k} - xbar);
    x{k+1} = A * x{k} + B * u(k);
end

cellx = cell2mat(x);

plot(u, 'k', 'LineWidth', 1);
xlabel('$k$', 'Interpreter','latex')
legend('$u(k) \ [N]$', 'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','FontSize',18);
figure; grid on;
yyaxis left
plot(cellx(1, :),  'LineWidth', 1);
ylim([-0.3, 1.5])
xlim([0, n_samples]);
yyaxis right
plot(cellx(2, :),  'LineWidth', 1);
ylim([-0.1, 0.5])
xlabel('$k$', 'Interpreter','latex')
legend('$x_1(k)  \ [m]$', '$x_2(k)  \ [m/s]$', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',18);

figure; grid on;
yyaxis left
plot(cellx(3, :),  'LineWidth', 1);
ylim([-0.2, 0.1])
xlim([0, n_samples]);
yyaxis right
plot(cellx(4, :),  'LineWidth', 1);
ylim([-0.2, 0.1])
xlabel('$k$', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',18);
legend('$x_3(k)  \ [rad]$', '$x_4(k)  \ [rad/s]$', 'Interpreter','latex')


