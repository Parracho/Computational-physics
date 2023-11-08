clc
clear all
close all



t0=0;
tf=100;
F0=1.0;
v0=2;
y0=0;
h=0.01;
epislon=0.3;

t=t0:h:tf;
N=length(t);
v=zeros(1,N);
y=zeros(1,N);

v(1)=v0;
y(1)=y0;



fy=@(V) (V);
fv=@(T,Y,V) (-epislon*(Y^2-1)*V-Y+F0*cos(1.7*T));

for i=1:N-1
    r1y=fy(v(i));
    r1v=fv(t(i),y(i),v(i));
    
    r2y=fy(v(i)+r1v*(h/3));
    r2v=fv(t(i)+(h/3),y(i)+r1y*(h/3),v(i)+r1v*(h/3));
    
    r3y=fy(v(i)-r1v*(h/3)+r2v*h);
    r3v=fv(t(i)+(2*h/3),y(i)-r1y*(h/3)+r2y*h,v(i)-r1v*(h/3)+r2v*h);
    
    r4y=fy(v(i)+r1v*h-r2v*h+r3v*h);
    r4v=fv(t(i)+h,y(i)+r1y*h-r2y*h+r3y*h,v(i)+r1v*h-r2v*h+r3v*h);
    
    y(i+1)=y(i)+(h/8)*(r1y+3*r2y+3*r3y+r4y);
    v(i+1)=v(i)+(h/8)*(r1v+3*r2v+3*r3v+r4v);
end

figure(1)
plot(t,v)
figure(2)
plot(y,v)