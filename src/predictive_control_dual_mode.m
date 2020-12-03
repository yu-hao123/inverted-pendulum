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

Ts = 0.1;
n_samples = 100;

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

Pf = S;
Af = A - B*K;
phi = B*(K*xbar + ubar);
Gamma = [eye(4); -K];
h = [zeros(4,1); K*xbar + ubar];

% Control System Restrictions
Sx = [1 0 0 0; -1 0 0 0; 0 0 1 0; 0 0 -1 0];
%Sx = [1 0; -1 0; 0 1;0 -1];
bx = [1; 0.1; pi/10; pi/10];
Su = [1; -1];
bu = [1; 1];

Spsi = blkdiag(Sx, Su);
bpsi = [bx; bu];

max_iter = 100;
tol = 0.01;
[Sf, bf] = determine_oinf(Af,phi,Gamma,h,Spsi,bpsi,max_iter,tol);

N = 30;

state_equation = @(t,x,u) ...
    [x(2); m*g*x(3)/M + u/M; x(4); (M+m)*g*x(3)/(M*L) + u/(M*L)];
x{1} = x0;
trec = 0; xrec = x{1}';
for k = 1:n_samples
    un = optimization_u(A, B, Q, R, Pf, N, Sx, bx, ...
    Su, bu, Sf, bf, x{k}, xbar, ubar);
    u(k) = un(1);
    [tout,xout] = ode45(@(t,x)state_equation(t,x,u(k)),[(k-1)*Ts k*Ts],x{k});
    trec = [trec;tout(2:end)];
    xrec = [xrec;xout(2:end,:)];
    x{k+1} = xrec(end,:)';
end

k = 0:n_samples;
figure;
stairs(k*Ts,[u u(end)], 'k', 'LineWidth', 1);
xlabel('$t \ [s]$', 'Interpreter', 'latex');
set(gca,'TickLabelInterpreter','latex','FontSize',18);
legend('$u(k) \ [N]$', 'Interpreter', 'latex')
ylim([-2, 2])
yline(1, 'r--', 'LineWidth', 1, 'HandleVisibility', 'off');
yline(-1, 'r--', 'LineWidth', 1, 'HandleVisibility', 'off');

figure, plot(trec,xrec(:,1), 'k', 'LineWidth', 1)
set(gca,'TickLabelInterpreter','latex','FontSize',18);
legend('$x_1(k)  \ [m]$', 'Interpreter','latex')
xlabel('$t \ [s]$', 'Interpreter', 'latex');
figure, plot(trec,xrec(:,2),  'k', 'LineWidth', 1)
set(gca,'TickLabelInterpreter','latex','FontSize',18);
legend('$x_2(k)  \ [m/s]$', 'Interpreter','latex')
xlabel('$t \ [s]$', 'Interpreter', 'latex');
ylim([-0.1, 0.4])
figure, plot(trec,xrec(:,3),  'k', 'LineWidth', 1)
set(gca,'TickLabelInterpreter','latex','FontSize',18);
legend('$x_3(k)  \ [rad]$', 'Interpreter','latex')
xlabel('$t \ [s]$', 'Interpreter', 'latex');
ylim([-0.1, 0.1])
figure, plot(trec,xrec(:,4),  'k', 'LineWidth', 1)
set(gca,'TickLabelInterpreter','latex','FontSize',18);
legend('$x_4(k)  \ [rad/s]$', 'Interpreter','latex')
xlabel('$t \ [s]$', 'Interpreter', 'latex');
ylim([-0.1, 0.1])