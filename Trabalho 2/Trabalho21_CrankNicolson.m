clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Crank Nicolson %%%%%%%%%%%%%%%%%%%%%

K=1;
m=1; 
w=sqrt(K/m);
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
    x(i+1)=(x(i)*(1-(w*h/2)^2)+h*v(i))/(1+(w*h/2)^2);
    v(i+1)=v(i)-(K/m)*x(i)*h;
end
E=(1/2)*(K*x.^2+m*v.^2);

k=1;
for j=2:N-1
    if x(j)>x(j-1) && x(j)>x(j+1)
        max_t(k)=t(j);
        max_x(k)=x(j);
        k=k+1;
    end
end

A_crank=mean(max_x)

T=0;
for i=1:k-2
    T=T+(max_t(i+1)-max_t(i));
end

T_crack = T/(double(k))


figure(1)
plot(t,x,"-r") %TEMPO-X POSICAO-Y
xlabel("Tempo(s)")
ylabel("Posição(m)")
title("CRANK-NICOLSON Method - Mola")
figure(2)
plot(t,v,"-r") %TEMPO-X VELOCIDADE-Y
xlabel("Tempo(s)")
ylabel("Velocidade(m/s)")
title("CRANK-NICOLSON Method - Mola")
figure(3)
plot(t,E,"-r") %TEMPO-X VELOCIDADE-Y
xlabel("Tempo(s)")
ylabel("Total Energy")
title("CRANK-NICOLSON Method - Mola")