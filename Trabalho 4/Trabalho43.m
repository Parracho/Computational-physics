clc
clear all
close all

epsilon=1;
t0=0;
tf=100;

tspan=[t0 tf];
y0=[1 2 3];
v0=[3 2 1];


n=length(y0);

for i=1:n
    options=odeset("RelTol",3e-14,"AbsTol",[1e-13 1e-13]);
    [t, solucao]=ode45(@f,tspan,[y0(i) v0(i)],options,epsilon);
    y=solucao(:,1);
    v=solucao(:,2);
    
    ts=round(length(t)/2);
    t=t(ts:end);
    
    plot(y,v)
    hold on
end

function derivadas = f(t,solucao,epsilon)
    f1=@(Y,V) (V);
    f2=@(Y,V) (-epsilon*(Y^2-1)*V-Y);
    derivadas=zeros(2,1);
    derivadas(1)=f1(solucao(1),solucao(2));
    derivadas(2)=f2(solucao(1),solucao(2));
end