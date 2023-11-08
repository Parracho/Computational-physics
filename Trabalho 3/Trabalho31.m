clc
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3.1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=1;
k=16;
w=sqrt(k/m);
v0=0;
x0=1;
t0=0;
h=0.10;
tf=12;
t=t0:h:tf;
N=length(t);
x=zeros(1,N);
v=zeros(1,N);
x(1)=x0;
v(1)=v0;
fx = @(V) (V);
fv = @(X) (- k/m*X);

for i=1:N-1
    r1x=fx(v(i));
    r1v=fv(x(i));
    r2x=fx(v(i)+r1v*h/2);
    r2v=fv(x(i)+r1x*h/2);
%     r1x=v(i);
%     r1v=-w^2*x(i);
%     r2x=v(i)+r1v*h/2;
%     r2v=-w^2*(x(i)+r1x*h/2);
    v(i+1)=v(i)+h*r2v;
    x(i+1)=x(i)+h*r2x;
end
E=(1/2)*(k*(x.^2)+m*(v.^2));

figure(1)
plot(t,x)
figure(2)
plot(t,v)
figure(3)
plot(t,E)















































