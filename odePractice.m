clear;clc; close all;

tspan = [0 5];
y0 = 0;
[t,y] = ode45(@(t,y) 2*t, tspan, y0);

figure
plot(t,y,'-o')


%%
[t,y] = ode45(@vdp1,[0 20],[2; 0]);

figure
plot(t,y(:,1),'-o',t,y(:,2),'-o')
title('Solution of van der Pol Equation (\mu = 1) with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2')

%%

addpath("func\");

A = 1;
B = 2;
tspan = [0 5];
y0 = [0 0.01];
[t,y] = ode45(@(t,y) odefcn(t,y,A,B), tspan, y0);

figure
plot(t,y(:,1),'-o',t,y(:,2),'-.')


%%
muk = 0; sigk = 1;

tspan = [0 20]; ic = 0.001;
opts = odeset(RelTol=1e-2,AbsTol=1e-4);
[t,x] = ode45(@(t,x) andersson2(x,muk,sigk), tspan,ic,opts);

figure;
plot(t,x);
grid on

%% Deprecated from modelThings.m

F = @(t,x) -mu_k*x;         % Eqn 1 in Andersson (2021); main equation in Campbell (1995)
C1 = modelEquation(F,t0,h,tF,x0);
