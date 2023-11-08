clc
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Analitical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=1;
m=1;
w=sqrt(k/m);        
x0=1;
v0=0;
t0=0;
tf=30;
h=0.1;
t=t0:h:tf; %VETOR TEMPO
N=length(t);
v=zeros(1,N);
x=zeros(1,N);
E=zeros(1,N-1);
v(1)=v0;
x(1)=x0;

for i=1:N-1
    v(i+1)=-x0*w*sin(w*t(i));
    x(i+1)=x0*cos(w*t(i));
end
E=(1/2)*(k*x.^2+m*v.^2);

figure(1)
plot(t,x,"-r") %TEMPO-X POSICAO-Y
xlabel("Tempo(s)")
ylabel("Posição(m)")
title("Analitical Method - Mola")
figure(2)
plot(t,v,"-r") %TEMPO-X VELOCIDADE-Y
xlabel("Tempo(s)")
ylabel("Velocidade(m/s)")
title("Analitical Method - Mola")
figure(3)
plot(t,E,"-r") %TEMPO-X VELOCIDADE-Y
xlabel("Tempo(s)")
ylabel("Total Energy")
title("Analitical Method - Mola")


