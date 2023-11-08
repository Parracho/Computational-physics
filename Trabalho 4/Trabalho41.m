clc
clear all
close all



t0=0;
tf=100;
h=0.01;
epsilon=0.7;
y0=0.2;
v0=0.7;
t=t0:h:tf;
N=length(t);
v=zeros(1,N);
y=zeros(1,N);
v(1)=v0;
y(1)=y0;

fy= @(V) (V);
fv= @(Y,V) (-epsilon*(Y^2-1)*V-Y);

for i=1:N-1
    r1y=fy(v(i));
    r1v=fv(y(i),v(i));
    
  
    r2y=fy(v(i)+r1v*h);
    r2v=fv(y(i)+r1y*h,v(i)+r1v*h);
    
    r3y=fy(v(i)+(r1v+r2v)*(h/4));
    r3v=fv(y(i)+(r1y+r2y)*(h/4),v(i)+(r1v+r2v)*(h/4));
    
    v(i+1)=v(i)+(h/6)*(r1v+r2v+5*r3v);
    y(i+1)=y(i)+(h/6)*(r1y+r2y+5*r3y);


end

ts=round(N/2);
t=t(ts:end);
v=v(ts:end);
y=y(ts:end);

plot(t,y)