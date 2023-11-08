clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Euler-Cromer  %%%%%%%%%%%%%%%%%%%%%%%%
h=0.001;
t=0:h:1;
N=length(t);
r=zeros(2,N);
rs=zeros(N);
v=zeros(2,N);
ang=zeros(N);
a=zeros(2,N);

r(1,1)=0.47;  %%%   y=0
rs(1)=0.47;
v(2,1)=8.2;   %%%   vx=0
ang(1)=0;



for i=1:N-1
    ac=-(norm(v(:,i)))^2/rs(i);
    a(1,i)=ac*cos(ang(i));
    a(2,i)=ac*sin(ang(i));
    v(:,i+1)=v(:,i)+a(:,i)*h;
    r(:,i+1)=r(:,i)+v(:,i+1)*h; 
    rs(i+1)=norm(r(:,i+1));
    ang(i+1)=atan(r(2,i+1)/r(1,i+1));
end

plot(r(1,:),r(2,:),"-r")
axis([-0.5 0.5 -0.5 0.5])
set(gca,'PlotBoxAspectRatio',[1 1 1])

