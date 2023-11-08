clc
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  1.2    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v0=50;
theta_d=37;
theta=theta_d*2*pi/360;  %degrees para radians
x0=0;
y0=0; %% alinea d)
%y0=55;
g=[0;-9.8];
m=1;
t0=0;
tf=10;
h=0.01;
t=t0:h:tf;
N=length(t);
r=zeros(2,N);
v=zeros(2,N);
r(1,1)=x0;
r(2,1)=y0;
vx0=v0*cos(theta);
vy0=v0*sin(theta);
v(1,1)=vx0;
v(2,1)=vy0;

for i=1:N-1
    v(:,i+1)=v(:,i)+g*h;  %vx and vy
    r(:,i+1)=r(:,i)+v(:,i+1)*h;    %x and y
    if r(2,i+1)<0   %%y=0
        break
    end
end


R_obtido=interp1(r(2,i:i+1),r(1,i:i+1),0,"linear");

R_teorico=-v0^2*sin(2*theta)/g(2)
R_obtido   %%y0=0 
figure(1)
plot(r(1,1:i),r(2,1:i),".r")
title("2D")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  1.3_centrifugal %%%%%%%%%%%%%%%%%%%%%%%%%%%



alpha=40.6;
R=6.37*10^6;
w=7.292*10^-5;
H=200;
m=1;
t0=0;
tf=10;
h=0.01;
t=t0:h:tf;
N=length(t);
v=zeros(3,N);
r=zeros(3,N);
g=[0;0;-9.8];
r(3,1)=H;
w=[-w*cos(alpha);0;w*sin(alpha)];
R=[0 0 R];
a_centrifugal=cross(w,cross(w,R));
a=g-a_centrifugal.';

for i=1:N-1
    v(:,i+1)=v(:,i)+a*h;
    r(:,i+1)=r(:,i)+v(:,i+1)*h;
    if r(3,i+1)<0
        break
    end
end

x_cent=interp1(r(3,i:i+1),r(1,i:i+1),0,"linear")
y_cent=interp1(r(3,i:i+1),r(2,i:i+1),0,"linear")

figure(2)
plot3(r(1,1:i),r(2,1:i),r(3,1:i),".r")
title("3D-centrifugal")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%    1.3_coriolis %%%%%%%%%%%%%%%%%%%%%%%%%%%


alpha=40.6;
R=6.37*10^6;
w=7.292*10^-5;
H=200;
m=1;
t0=0;
tf=10;
h=0.01;
t=t0:h:tf;
N=length(t);
v=zeros(3,N);
r=zeros(3,N);
a=zeros(3,N);
g=[0;0;-9.8];
r(3,1)=H;
a(:,1)=g;
w=[-w*cos(alpha);0;w*sin(alpha)];
R=[0 0 R];

for i=1:N-1
    v(:,i+1)=v(:,i)+a(:,i)*h;
    a(:,i+1)=g-2*cross(w,v(:,i+1));
    r(:,i+1)=r(:,i)+v(:,i+1)*h;
    if r(3,i+1)<0
        break
    end
end

x_cori=interp1(r(3,i:i+1),r(1,i:i+1),0,"linear")
y_cori=interp1(r(3,i:i+1),r(2,i:i+1),0,"linear")

figure(3)
plot3(r(1,1:i),r(2,1:i),r(3,1:i),".r")
title("3D-coriolis")










































