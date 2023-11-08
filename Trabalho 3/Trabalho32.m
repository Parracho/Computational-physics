clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% a,b,c) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m=1;
% k=16;
% w=sqrt(k/m);
% v0=0;
% x0=1;
% t0=0;
% h=0.3;
% tf=12;
% t=t0:h:tf;
% N=length(t);
% x=zeros(1,N);
% v=zeros(1,N);
% xa=zeros(1,N);
% va=zeros(1,N);
% x(1)=x0;
% v(1)=v0;
% fx = @(V)V;
% fv = @(X) - k/m*X;
% 
% for i=1:N %para terem o msm tamanho dps "elimino" no plot
%     r1x=fx(v(i));
%     r1v=fv(x(i));
%     r2x=fx(v(i)+r1v*h/2);
%     r2v=fv(x(i)+r1x*h/2);
%     r3x=fx(v(i)+r2v*h/2);
%     r3v=fv(x(i)+r2x*h/2);
%     r4x=fx(v(i)+r3v*h);
%     r4v=fv(x(i)+r3x*h);
%     v(i+1)=v(i)+(h/6)*(r1v+2*r2v+2*r3v+r4v);
%     x(i+1)=x(i)+(h/6)*(r1x+2*r2x+2*r3x+r4x);
%     xa(i)=x0*cos(w*t(i));
%     va(i)=-w*x0*sin(w*t(i));
% end
% E=(1/2)*(k*(x.^2)+m*(v.^2));
% Ea=(1/2)*(k*(x0^2)+m*(v0^2));
% 
% figure(1)
% plot(t,x(1:end-1),"-b")
% hold on
% plot(t,xa,"-r")
% figure(2)
% plot(t,v(1:end-1),"-b")
% hold on
% plot(t,va,"-r")
% figure(3)
% plot(t,E)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% d) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m=1;
% k=16;
% w=sqrt(k/m);
% v0=0;
% x0=1;
% t0=0;
% tf=10;
% errog=[];
% fx = @(V)V;
% fv = @(X) - k/m*X;
% h=[0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2];
% for j=1:length(h)
%     
%     t=t0:h(j):tf;
%     N=length(t);
%     x=zeros(1,N);
%     v=zeros(1,N);
%     xa=zeros(1,N);
%     va=zeros(1,N);
%     x(1)=x0;
%     v(1)=v0;
%     
% 
%     for i=1:N 
%         r1x=fx(v(i));
%         r1v=fv(x(i));
%         r2x=fx(v(i)+r1v*h(j)/2);
%         r2v=fv(x(i)+r1x*h(j)/2);
%         r3x=fx(v(i)+r2v*h(j)/2);
%         r3v=fv(x(i)+r2x*h(j)/2);
%         r4x=fx(v(i)+r3v*h(j));
%         r4v=fv(x(i)+r3x*h(j));
%         v(i+1)=v(i)+(h(j)/6)*(r1v+2*r2v+2*r3v+r4v);
%         x(i+1)=x(i)+(h(j)/6)*(r1x+2*r2x+2*r3x+r4x);
%         xa(i)=x0*cos(w*t(i));
%         va(i)=-w*x0*sin(w*t(i));
%     end
%     errog(j)=log10(abs(x(end-1)-xa(end)));
% end
% coeff=polyfit(log10(h),errog,1)
% y=coeff(1)*log10(h)+coeff(2);
% 
% figure(1)
% plot(log10(h),errog,"-r")
% hold on
% plot(log10(h),y,"-b")
% 














































