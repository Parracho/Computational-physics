clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    a)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dy/dt=v
% dv/dt=(miu/m)*sin(v)*v+(F0/m)*cos(w0*t)-(K/m)*(y+a*y^3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    b)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=1;
K=1;
a=2;
miu=2;
t0=0;
tf=150;
tspan=[t0 tf];
y0=1.5;
v0=0;
F0=[0 0.5 1.5];
n=length(F0);
w0=1.5;


for i=1:n
    options=odeset("RelTol",3e-14,"AbsTol",[1e-13 1e-13]);
    [t, solucao]=ode45(@(t,y,v) f(t,y,F0(i),miu,a,K,m,w0),tspan,[y0 v0],options);
    y=solucao(:,1);
    v=solucao(:,2);
    
    ts=round(length(t)/2);
    t=t(ts:end);
    y=y(ts:end);
    v=v(ts:end);
    
    figure(i)
    plot(t,y)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    d)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=1;
K=1;
a=2;
miu=2;
y0=1.5;
v0=0;
t0=0;
tf=150;
h=0.001;
t=t0:h:tf;
N=length(t);
y=zeros(1,N);
v=zeros(1,N);
v(1)=v0;
y(1)=y0;
F0=0;
w0=1.5;

fy=@(V)(V);
fv=@(T,Y,V)((miu/m)*sin(V)*V+(F0/m)*cos(w0*T)-(K/m)*(Y+a*Y^3));

for i=1:N-1
    r1y=fy(v(i));
    r1v=fv(t(i),y(i),v(i));

    r2y=fy(v(i)+r1v*h/3);
    r2v=fv(t(i)+h/3,y(i)+r1y*h/3,v(i)+r1v*h/3);

    r3y=fy(v(i)-r1v*h/3+r2v*h);
    r3v=fv(t(i)+2*h/3,y(i)-r1y*h/3+r2y*h,v(i)-r1v*h/3+r2v*h);

    r4y=fy(v(i)+r1v*h-r2v*h+r3v*h);
    r4v=fv(t(i)+h,y(i)+r1y*h-r2y*h+r3y*h,v(i)+r1v*h-r2v*h+r3v*h);

    y(i+1)=y(i)+(h/8)*(r1y+3*r2y+3*r3y+r4y);
    v(i+1)=v(i)+(h/8)*(r1v+3*r2v+3*r3v+r4v);
end

figure(4)
plot(t,y)


function derivadas = f(t,solucao,F0,miu,a,K,m,w0)
    f1=@(V) (V);
    f2=@(T,Y,V) ((miu/m)*sin(V)*V+(F0/m)*cos(w0*T)-(K/m)*(Y+a*Y^3));
    derivadas=zeros(2,1);
    derivadas(1)=f1(solucao(2));
    derivadas(2)=f2(t,solucao(1),solucao(2));
end
















