clc
clear all
close all


t0=0;
tf=10;
x0=1;
v0=0;
tspan=[t0 tf];
y0=[x0 v0];
m=1;
k=16;

options=odeset('RelTol',3e-14,'AbsTol',[1e-13 1e-13]);
[t,y]=ode45(@odefun, tspan, y0, options, m, k);
x=y(:,1);
v=y(:,2);
E=(1/2)*(k*(x.^2)+m*(v.^2));
figure(1)
plot(t,y(:,1),".r",t,y(:,2),".b")
figure(2)
plot(t,E)
figure(3)
plot(x,v)
function dydt=odefun(t,y,m,k)
    dydt=zeros(2,1);
    dydt=[y(2);-(k/m)*y(1)];
end

