% script for the Runge-Kutta calculation of the dynamics
clc;
clear all;
global A1 A2 w1 w2 w0 kappa
r=1.0;
A1=0.2;
A2=r*A1;
w1=1.00;w2=1.005*w1;
w0=1.0;
kappa=0.0;
tspan=0:0.25:2000;
y0=[0;0;-1];
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[tt,yy]=ode45(@fundydt,tspan,y0,opts);
plot(tt,real(0.5*(1+yy(:,3)))),hold on
% [tt2,yy2]=ode45(@fundydtrwa,tspan,y0,opts);
% plot(tt2,0.5*(1+yy2(:,3))),hold on