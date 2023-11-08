clc;
clear all;
close all;
h=logspace(-4,-1,10);
errog=[];
for j=1:length(h)
z0=6;
v0=0;
t0=0;
g=9.8;
%h=0.2;
t=t0:h:1.5;
N=length(t);
v=zeros(1,N);
z=zeros(1,N);
z(1)=z0;
v(1)=v0;
m=0.15;

for i=1:N
    if z(i)>=0
        v(i+1)=v(i)-g*h;
        z(i+1)=z(i)+v(i)*h;
    end
end
zmax=max(z);
errog(j)=abs(z0-zmax);
end
plot(ln(h),ln(errog))