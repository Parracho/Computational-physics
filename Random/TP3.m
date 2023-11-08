% %TP3
% 
 close all
 clear all
 clc
% 
% %2.1
% 
% t0 = 0;
% tf = 500;
% h = 0.1; %dt
% 
% m = 1; %massa
% K = 1; %constante da mola
% x0 = 1; %deslocamento inicial
% 
% t = (t0:h:tf);
% k = length(t);
% 
% 
% vx = (1:k);
% vx(1)= 0;
% x = (1:k);
% x(1) = x0;
    

% %Euler
% for j = 1:(k-1)
% 
%     vx(j+1) = vx(j) - (K/m) * x(j) * h;
%     x(j+1) = x(j) + vx(j) * h; 
%    
% end


% %Euler-Cromer
% for j = 1:(k-1)
%     
%     vx(j+1) = vx(j) - (K/m) * x(j) * h;
%     x(j+1) = x(j) + vx(j+1) * h;
% 
% end


% % Euler Implicito
% for j = 1:(k-1)
%     
%     x(j+1) = (x(j) + vx(j) * h) / (1 + (K/m) * h^2); 
%     vx(j+1) = vx(j) - (K/m) * x(j+1) * h; 
%     
% end

% %Cranck-Nicolson
% for j = 1:(k-1)
%     
%     x(j+1) = (x(j) + h * vx(j) - (h^2 /4)* (K/m) * x(j))/(1+(K/m)*(h^2 / 4));
%     vx(j+1) = vx(j) - (K/m) * (x(j) + x(j+1)) * (h/2); 
%     
% end
% 
% l = 1; %L
% 
% for i = 2:k-1
%     
%     if x(i+1)-x(i)<=0 && x(i)-x(i-1)>=0
%         aux = lagr(t(i-1:i+1),x(i-1:i+1));
%         xmax(l) = aux(2);
%         tmax(l) = aux(1); 
%         l = l+1;
%     end
%     
% end
% 
% tmax_polyfit = polyfit(1:l-1,tmax,1)
% 
% periodo_average = tmax_polyfit(1)
% 
% E = (1/2) * K*(x.^2) + (1/2) * m * (vx.^2);
% 
% amp_average = mean(xmax);
% p_average = mean(tmax);

% plot(t,x,'b')
% figure(2)
% plot(t,vx,'r')
% figure(3)
% plot(t,E,'k')
% figure(4)
% plot(x,vx,'b')



%2.2

% h = 0.0001;
% t0 =0;
% tf= 10; %anos
% t=(t0:h:tf);
% k = length(t);
% 
% x0 = 0.47; %AU - unidades astronomicas
% y0 = 0;
% vx0 = 0;
% vy0 = 8.2; %AU/ano
% 
% angulo = 0; %graus
% 
% %v =zeros(2,k);
% v(1,1) = vx0;
% v(2,1) = vy0;
% 
% %r = zeros(2,k);
% r(1,1) = x0;
% r(2,1) = y0;
% 
% for i = 1:k-1
%     
%     raio = sqrt(r(1,i)^2 + r(2,i)^2); %raio em cada ponto
%     angulo = mod(atan2(r(2,i),r(1,i)),2*pi);
%     
%     v(1,i+1) = v(1,i) - ((4*(pi^2) / raio^3) * (raio * cos(angulo))) * h;
%     v(2,i+1) = v(2,i) - ((4*(pi^2) / raio^3) * (raio * sin(angulo))) * h; %velocidades
%     
%     r(1,i+1) = r(1,i) + v(1,i+1) * h;
%     r(2,i+1) = r(2,i) + v(2,i+1) * h; %posiçoes
%     
%     if r(1,i)
%         break
%     end
%     
% end
% 
% plot(r(1,:),r(2,:));
% hold on
% 
% axis([-0.5 0.5 -0.5 0.5])
% set(gca,'PlotBoxAspectRatio',[1 1 1])


%%%%%%2.3

h = 0.01;
t0 =0;
tf= 10; 
t=(t0:h:tf);
k = length(t);

x0 = 1;
vx0 = 1;
alfa = -0.1;
K = 1;
m = 1;

x = zeros(1,k);
x(1) = x0;
vx = zeros(1,k);
vx(1)=vx0;

V = zeros(1,k); %energia potencial  

    for i = 1:k-1

       vx(i+1) = vx(i) - (K/m) *(x(i)+(2*alfa * (x(i)^3))) * h;
       x(i+1) = x(i) + vx(i+1) * h;

       V(i) = (1/2) * K * ((x(i)^2)+ alfa * (x(i)^4));

    end

figure(1)
subplot(1,3,1)
plot(t,x,'b')
subplot(1,3,2)
plot(t,vx,'r')
subplot(1,3,3)
plot(t,V,'k')

