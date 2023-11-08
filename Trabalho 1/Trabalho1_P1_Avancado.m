clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%  1.1.a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% g=9.8;
% vlim=6.8;
% t0=0;
% tf=2;
% v0=0;
% h=0.01;
% t=t0:h:tf;
% N=length(t);
% vz=[v0];
% va=[v0];
% f=[];
% for j=1:N-1
%     va(j+1)=-vlim*tanh(g*t(j+1)/vlim);
%     f(j)=-g*(1-(tanh(g*t(j)/vlim))^2); %% f(t)=va'(t)
%     vz(j+1)=vz(j)+h*f(j);
% end
% figure
% plot(t,va,"-r")
% hold on
% plot(t,vz,"-b")

%%%%%%%%%%%%%%%%%%%%%%%%%%  1.1.b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% g=9.8;
% vlim=6.8;
% v0=0;
% h=0.01:0.01:0.1;
% erro=zeros(1,length(h));
% t=0.5;
% for i=1:length(h)
%     va=-vlim*tanh(g*t/vlim);
%     f=-g*(1-(tanh(g*t/vlim))^2); %% f(t)=va'(t)
%     vz=v0+h(i)*f;
%     erro(i)=abs(vz-va); 
% end
% figure(1)
% plot(log10(h),log10(erro),"rx:")
% xlabel("log Precisao")
% ylabel("log Erro Global")
% figure(2)
% plot(h,erro,"rx:")

%%%%%%%%%%%%%%%%%%%%%%%%%%  1.1.c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


g=9.8;
vlim=6.8;
t0=0;
tf=100;
v0=16;
z0=1;
h=0.1;
t=t0:h:tf;
N=length(t);
vz=zeros(1,N-1);
zz=zeros(1,N-1);
vz(1)=v0;
zz(1)=z0;


for j=1:N-1
    vz(j+1)=vlim^2*(1/(g*t(j+1)));
    zz(j+1)=zz(j)+h*vz(j+1);
end


for i=1:N-1  % determinar o zero e fazer a interpolação
    if zz(i)*zz(i+1)<0 
        t_gz=interp1(zz(i:i+1),t(i:i+1),0,"linear");
        break
    end
end

figure

plot(t(1:i),vz(1:i),"-r")




%%%%%%%%%%%%%%%%%%%%%%%%%%  1.2.  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%  1.2.a  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



























