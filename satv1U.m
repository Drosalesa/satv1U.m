%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Balance térmico en un CubeSat
%
% 08/04/21

%% Datos

l = 0.1;
a = 0.1;
al = 0.1;
Is = 1367;
alpha_u = 0.5;
alpha_c = 0.91;
eff_c = 0.28;
fc = 0.8;
Re = 6378e3;
R = 7110e3;
%% Paso y vector v de angulos

p = 1000;
v(1:p+1,1)=0;
k=0;
for i = 1:p+1
    v(i)= k;
    k=k+2*pi()/p;
end
%% Areas

A(1:6,1)=0;
A(1)=a*al;
A(2:5) = al*l;
A(6)=a*al;
%% Alpha promedio

alpha_eff = alpha_c-eff_c;
alpha(1:6,1)=0;
alpha(1)= alpha_u;
alpha(2:5)= fc*alpha_eff+(1-fc)*alpha_u;
alpha(6)=alpha_u;
%% Radiacion directa del sol

Qs(1:p+1,1)=0;
Qsi(1:6,1)=0;
eps = asin(Re/R);
B(1:6,1)=0;
for j=1:p+1
    for i = 1:6
    B(1)=-cos(pi()/2-v(j));
    if B(1)<=0
        B(1)=0;
    end
    B(2)=-cos(v(j));
    if B(2)<=0
        B(2)=0;
    end
    B(3)=0;
    B(4)=cos(v(j));
    if B(4)<=0
        B(4)=0;
    end
    B(5)=0;
    B(6)=cos(pi()/2-v(j));
    if B(6)<=0
        B(6)=0;
    end
    Qsi(i)=alpha(i)*A(i)*B(i)*Is;
    if v(j)>(pi()-eps) && v(j)<(pi()+eps)
        Qsi(i)=0;
    end
    end
    Qs(j)=sum(Qsi);
end
%% Albedo
gamma=0.273;
h = R/Re;
k=0;
the_m=acos(Re/R);
Qalb1(1:p+1,1)=0;
Qalb2(1:p+1,1)=0;
Qalb3(1:p+1,1)=0;
v(1:p+1,1)=0;
for i=1:p+1
    v(i)=k;
    k=k+2*pi()/p;
end

s1 = (pi()/2-the_m)/(2*pi()/(p+1));
s1 = round(s1); % pi()/2 - the_m
s2 = (pi()/2)/(2*pi/(p+1)); 
s2 = round(s2); % % pi()/2
s3  = (pi()/2+the_m)/(2*pi/(p+1));
s3 = round (s3); % pi()/2 + the_m
s4 = (3*pi()/2-the_m)/(2*pi()/(p+1));
s4 = round (s4); % 3*pi()/2 - the_m
s5 = (3*pi()/2)/(2*pi()/(p-1));
s5 = round(s5); % 3*pi()/2
s6 = (3*pi()/2+the_m)/(2*pi()/(p+1));
s6 = round(s6); % 3*pi()/2 + the_m

for i=1:p+1
    F1Ea = @(x,y) (1-h.*cos(x)).*cos(y).*sin(x)./(1+h.^2-2.*h.*cos(x)).^2.*(sin(x).*cos(y).*sin(v(i))+cos(x).*cos(v(i))).*sin(x);
    F2Ea = @(x,y) (1-h.*cos(x)).*(cos(x)-h).*sin(x)./(1+h.^2-2.*h.*cos(x)).^2.*(sin(x).*cos(y).*sin(v(i))+cos(x).*cos(v(i)));
    F3Ea = @(x,y) (1-h.*cos(x)).*sin(y).*sin(x)./(1+h.^2-2.*h.*cos(x)).^2.*(sin(x).*cos(y).*sin(v(i))+cos(x).*cos(v(i))).*sin(x);
    phi_a = @(x) acos(-cot(v(i))*cot(x));
    phi_b = @(x) 2*pi()-acos(-cot(v(i))*cot(x));
    if i <= s1
        Qalb1(i) = 1/pi()*quad2d(F1Ea,0,the_m,pi()/2,3*pi()/2);
        Qalb2(i) = 1/pi()*quad2d(F2Ea,0,the_m,0,2*pi());
        Qalb3(i) = 1/pi()*quad2d(F3Ea,0,the_m,pi(),2*pi());
    end
    if i > s1 && i<=s2
        Qalb1(i) = 1/pi()*(quad2d(F1Ea,0,pi()/2-v(i),pi()/2,3*pi()/2)+quad2d(F1Ea,pi()/2-v(i),the_m,pi()/2,phi_a)+quad2d(F1Ea,pi()/2-v(i),the_m,phi_b,3*pi()/2));
        Qalb2(i) = 1/pi()*(quad2d(F2Ea,0,pi()/2-v(i),0,2*pi())+quad2d(F2Ea,pi()/2-v(i),the_m,0,phi_a)+quad2d(F2Ea,pi()/2-v(i),the_m,phi_b,2*pi()));
        Qalb3(i) = 1/pi()*(quad2d(F3Ea,0,pi()/2-v(i),pi(),2*pi())+quad2d(F3Ea,pi()/2-v(i),the_m,phi_b,2*pi()));
    end
    if i > s2 && i<=s3
        Qalb1(i) = 0;
        Qalb2(i) = 1/pi()*(quad2d(F2Ea,v(i)-pi()/2,the_m,0,phi_a)+quad2d(F2Ea,v(i)-pi()/2,the_m,phi_b,2*pi()));
        Qalb3(i) = 1/pi()*(quad2d(F3Ea,v(i)-pi()/2,the_m,phi_b,2*pi()));
    end
    if i > s3 && i<=s4
        Qalb1(i) = 0;
        Qalb2(i) = 0;
        Qalb3(i) = 0;
    end
    if i > s4 && i<=s5
        Qalb1(i) = 1/pi()*(quad2d(F1Ea,3*pi()/2-v(i),the_m,phi_a,phi_b));
        Qalb2(i) = 1/pi()*(quad2d(F2Ea,3*pi()/2-v(i),the_m,phi_a,phi_b));
        Qalb3(i) = 1/pi()*(quad2d(F3Ea,3*pi()/2-v(i),the_m,pi(),phi_b));
    end
    if i > s5 && i<=s6
        Qalb1(i) = 1/pi()*(quad2d(F1Ea,0,the_m,pi()/2,3*pi()/2));
        Qalb2(i) = 1/pi()*(quad2d(F2Ea,0,v(i)-3*pi()/2,0,2*pi())+quad2d(F2Ea,v(i)-3*pi()/2,the_m,phi_a,phi_b));
        Qalb3(i) = 1/pi()*(quad2d(F3Ea,0,v(i)-3*pi()/2,pi(),2*pi())+quad2d(F3Ea,v(i)-3*pi()/2,the_m,pi(),phi_b));
    end
    if i > s6
        Qalb1(i) = 1/pi()*(quad2d(F1Ea,0,the_m,pi()/2,3*pi()/2));
        Qalb2(i) = 1/pi()*(quad2d(F2Ea,0,the_m,0,2*pi()));
        Qalb3(i) = 1/pi()*(quad2d(F3Ea,0,the_m,pi(),2*pi()));
    end
end
Qalb=(alpha(2).*Qalb2.*A(2)+2.*alpha(3).*Qalb3.*A(3)+2.*alpha(1).*Qalb3.*A(1)).*gamma.*Is;

v_ang(1:p+1,1)=0;
for i=1:p+1
    v_ang(i)=v(i)*180/pi();
end
%{
figure(1)
plot(v_ang,Qalb1);
grid on;
xticks([45 90 135 180 225 270 315 360])
figure(2)
plot(v_ang,Qalb2); grid
xticks([45 90 135 180 225 270 315 360])
figure(3)
plot(v_ang,Qalb3); grid
xticks([45 90 135 180 225 270 315 360])
pause();
%}
%% Radiacion infrarroja
Ie = 213;
FiE(1:6,1)=0;
em_u=0.05;
em_c=0.89;
em_av=fc*em_c+(1-fc)*em_u;
F2 = @(x,y) (1-h.*cos(x)).*(cos(x)-h)./(1+h.^2-2.*h.*cos(x)).^2.*sin(x);
F6 = @(x,y) (h.*cos(x)-1).*sin(x).*cos(y)./(1+h.^2-2.*h.*cos(x)).^2.*sin(x);
for i=1:6
    FiE(i)=1/pi()*quad2d(F6,0,the_m,-pi()/2,pi()/2);
    if i==2
        FiE(i)=1/pi()*quad2d(F2,0,the_m,0,2*pi());
    end
    if i==4
        FiE(i)=0.0;
    end
end
Em(1:6,1)=0;
for i=1:6
    Em(i)=em_u;
    if i>1 && i<6
        Em(i)=em_av;
    end
end
QEi=alpha.*A.*FiE.*Ie;
QE(1:p+1,1)=sum(QEi);
%% Calor interno generado
QgenM = 0;
Qgenm = 0;
Qint(1:p+1,1)=Qgenm;
for i=1:p+1
    if v(i)>(pi()-eps) && v(i)<(pi()+eps)
        Qint(i)=QgenM;
    end
end

%% Velocidad
mu = 3.986004415e5;
w =(mu/((R/1000)^3))^(1/2);
P = 2*pi()/w/60;
dt = (2*pi/p)/w;

%% Solucion
Qtot=Qalb+QE+Qs+Qint;
C=336.34; %[J/°C]
sigma = 5.670373e-8;
A_eff=sum(A.*Em);

Ti=273;
Tn1(1:p+1,1)=0;
Tn1(1)=Ti;
for i=2:p+1
    Tn1(i)=Tn1(i-1)+dt/C*(Qtot(i-1)-sigma*A_eff*(Tn1(i-1))^4); 
end

%{
Tn2(1:p+1,1)=0;
for i=1:p+1
    Tn2(i)=(Qtot(i)/(2*sigma*A_eff))^(1/4)-273;
end
%}
%% Grafica

v_ang(1:p+1,1)=0;
for i=1:p+1
    v_ang(i)=v(i)*180/pi();
end
T_c(1:p+1,1)=0;
for i=1:p+1
    T_c(i)=Tn1(i)-273;
end
figure(1)
plot(v_ang,Qs);
hold on;
plot(v_ang,Qalb,'r --');
xticks([45 90 135 180 225 270 315 360])
axis([0 360 0 15]);
title('Radiacion directa del sol');
xlabel('v [°]');
ylabel('Q_{s}, Q_{alb} [W]');
grid on;
figure(2)
plot(v_ang,Qtot); grid
xticks([45 90 135 180 225 270 315 360])
title('Flujo de calor total');
xlabel('v [°]');
ylabel('Q_{tot} [W]');
figure(3)
plot(v_ang,T_c); grid
hold on;
%plot(v_ang,Tn2,'r');
xticks([45 90 135 180 225 270 315 360])
title('Temperatura de equilibrio');
xlabel('v [°]');
ylabel('T [°C]');
%axis([0 360 0 40]);

%% Temperatura de Estabilización
Tnf1(1:5*p+1)=0;
Qnf(1:5*p+1)=0;
time(1:5*p+1)=0;
for i=2:5*p+1
    time(i)=time(i-1)+dt/3600;
end
Tnf1(1)=Ti;
Tnf2(1)=Ti+15;
Tnf3(1)=Ti-15;
for i=1:p+1
    Qnf(i)=Qtot(i);
end
for i=p+2:2*p+1
    Qnf(i)=Qtot(i-p);
end
for i=2*p+2:3*p+1
    Qnf(i)=Qtot(i-2*p);
end
for i=3*p+2:4*p+1
    Qnf(i)=Qtot(i-3*p);
end
for i=4*p+2:5*p+1
    Qnf(i)=Qtot(i-4*p);
end
for i=2:5*p+1
    Tnf1(i)=Tnf1(i-1)+dt/C*(Qnf(i-1)-sigma*A_eff*(Tnf1(i-1))^4); 
end
for i=2:5*p+1
    Tnf2(i)=Tnf2(i-1)+dt/C*(Qnf(i-1)-sigma*A_eff*(Tnf2(i-1))^4); 
end
for i=2:5*p+1
    Tnf3(i)=Tnf3(i-1)+dt/C*(Qnf(i-1)-sigma*A_eff*(Tnf3(i-1))^4); 
end
T_cf1(1:5*p+1,1)=0;
T_cf2(1:5*p+1,1)=0;
T_cf3(1:5*p+1,1)=0;
for i=1:5*p+1
    T_cf1(i)=Tnf1(i)-273;
end
for i=1:5*p+1
    T_cf2(i)=Tnf2(i)-273;
end
for i=1:5*p+1
    T_cf3(i)=Tnf3(i)-273;
end
figure(4)
plot(time,T_cf1); grid
hold on;
plot(time,T_cf2,'r');
plot(time,T_cf3,'g');
title('Temperatura predicha del CubeSat');
xlabel('Tiempo [hr]');
ylabel('Temperatura [°C]');