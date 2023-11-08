clc;
clear all;
close all;
% %%###########################  Parte 1  ###############################
% 
% %%%%%%%%%%%%%%%%%%%%%%%%        A       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
z0=6; %POSICAO INICIAL
v0=0; %VELOCIDADE INICIAL
t0=0; %TEMPO INICIAL
g=9.8; %GRAVIDADE
h=0.2;
tf=1.5;
t=t0:h:tf; %VETOR TEMPO
N=length(t);
v=zeros(1,N); %VELOCIDADE
z=zeros(1,N); %POSICAO
z(1)=z0;
v(1)=v0;

for i=1:N-1
    v(i+1)=v(i)-g*h; %METODO DE EULER
    z(i+1)=z(i)+v(i+1)*h; %METODO DE EULER
end

for k=1:N-1  % determinar o zero e fazer a interpolação
    if z(k)*z(k+1)<0
        t_g=interp1(z(k:k+1),t(k:k+1),0,"linear");
        v_g=interp1(z(k:k+1),v(k:k+1),0,"linear");
        break
    end
end



figure(1)
plot(t(1:k),z(1:k),"-r") %TEMPO-X POSICAO-Y %% USEI t(1:j-1) para ignorar os termos que passam do Z=0
axis([0 t_g+0.2 -1 z0+1])
xlabel("Tempo(s)")
ylabel("Posição(m)")
title("Euler Method - Free Fall")
hold on
figure(2)
plot(t(1:k),v(1:k),"-r") %TEMPO-X VELOCIDADE-Y
axis([0 t_g+0.2 v_g-1 v0+1 ]) %Limites grafico
xlabel("Tempo(s)")
ylabel("Velocidade(m/s)")
title("Euler Method - Free Fall")

% 
% %#######################     ANALITICA        #########################
% 
z0=6; %POSICAO INICIAL
t0=0; %TEMPO INICIAL
g=9.8; %GRAVIDADE
h=0.2;
t=t0:h:1.5; %VETOR TEMPO
N=length(t);
z=zeros(1,N); %POSICAO
z(1)=z0;

for i=1:N-1
    z(i+1)=z(1)-(1/2)*g*t(i+1)^2;
end

for k=1:N-1  % determinar o zero e fazer a interpolação
    if z(k)*z(k+1)<0
        t_g=interp1(z(k:k+1),t(k:k+1),0,"linear");
        break
    end
end

figure(1)
plot(t(1:k),z(1:k),"-r") %TEMPO-X POSICAO-Y %% USEI t(1:j-1) para ignorar os termos que passam do Z=0
axis([0 t_g+0.2 -1 z0+1])
xlabel("Tempo(s)")
ylabel("Posição(m)")
title("Analitical Method - Free Fall")

%############################  Parte 2  ################################




%%###########################    METODO EULER      ###############################

k=2.5;
m=0.5;    
x0=10;
v0=0;
t0=0;
tf=10;
h=0.1;
t=t0:h:tf; %VETOR TEMPO
N=length(t);
v=zeros(1,N-1);
x=zeros(1,N-1);
v(1)=v0;
x(1)=x0;

for i=1:N-1
    v(i+1)=v(i)-(k/m)*x(i)*h;
    x(i+1)=x(i)+h*v(i);
end

figure(4)
plot(t,x,"-r") %TEMPO-X POSICAO-Y
xlabel("Tempo(s)")
ylabel("Posição(m)")
title("Euler Method - Mola")
figure(5)
plot(t,v,"-r") %TEMPO-X VELOCIDADE-Y
xlabel("Tempo(s)")
ylabel("Velocidade(m/s)")
title("Euler Method - Mola")


%%###########################   ANALITICA    #########################

k=2.5;
m=0.5;
w=sqrt(k/m);        
x0=10;
v0=0;
t0=0;
tf=10;
h=0.1;
t=t0:h:tf; %VETOR TEMPO
N=length(t);
v=[v0];
x=[x0];

for i=1:N-1
    v(i+1)=-x0*w*sin(w*t(i));
    x(i+1)=x0*cos(w*t(i));
end

figure(6)
plot(t,x,"-r") %TEMPO-X POSICAO-Y
xlabel("Tempo(s)")
ylabel("Posição(m)")
title("Analitical Method - Mola")
figure(7)
plot(t,v,"-r") %TEMPO-X VELOCIDADE-Y
xlabel("Tempo(s)")
ylabel("Velocidade(m/s)")
title("Analitical Method - Mola")


% %%%%%%%%%%%%%%%%%%%%%%%  Funcao Erro Global   %%%%%%%%%%%%%%%%%%%%%

k=2.5;
m=0.5;
w=sqrt(k/m);        
x0=10;
v0=0;
t0=0;
tf=10;
h=logspace(-4,-1,10);
erro_g=[];
for j=1:length(h)
    t=t0:h(j):tf; %VETOR TEMPO
    N=length(t);
    v=zeros(1,N-1);
    x=zeros(1,N-1);
    v(1)=v0;
    x(1)=x0;

    for i=1:N-1
        v(i+1)=v(i)-(k/m)*x(i)*h(j);
        x(i+1)=x(i)+h(j)*v(i);
    end
    xmax=max(x);
    erro_g(j)=abs(x0-xmax);
end
figure(8)
plot(h,erro_g,"-r")
xlabel("Precisão")
ylabel("Erro Global")
title("Precision vs Global Error")







































