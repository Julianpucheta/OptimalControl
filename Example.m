%Example extracted from: Sontag. Mathematical control theory 1998. Pag 104. http://www.sontaglab.org.
clc;clear all; clc;clear all;
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5;
h=0.0001;tiempo=(20/h);p_pp=0;tita_pp=0; t=0:h:tiempo*h;
%Condiciones iniciales
alfa(1)=.15; color='.g';
% alfa(1)=2.7; color='b';
alfa(1)=3.14; color='.r';
ref=0;
omega(1)=0; p_p(1)=0; u(1)=0; p(1)=0; i=1;indice=0;
%Versión linealizada en el equilibrio inestable. Sontag Pp 104.
% estado=[p(i); p_p(i); alfa(i); omega(i)]
Mat_A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(long*M) g*(m+M)/(long*M) 0]
Mat_B=[0; 1/M; 0; -1/(long*M)]
Mat_C=[1 0 0 0]; %La salida monovariable es posición y ángulo
Mat_M=[Mat_B Mat_A*Mat_B Mat_A^2*Mat_B Mat_A^3*Mat_B];%Matriz Controlabilidad
auto_val=eig(Mat_A);
c_ai=conv(conv(conv([1 -auto_val(1)],[1 -auto_val(2)]),[1 -auto_val(3)]),[1 -auto_val(4)]); %= poly(eig(Mat_A))
Mat_W=[;c_ai(4) c_ai(3) c_ai(2) 1;c_ai(3) c_ai(2) 1 0;c_ai(2) 1 0 0;1 0 0 0];
Mat_T=Mat_M*Mat_W;% T=Q=Co*Toeplitz
A_controlable=inv(Mat_T)*Mat_A*Mat_T; %Verificación de que T esté bien
a4=-A_controlable(4,1);%a4
a3=-A_controlable(4,2);
a2=-A_controlable(4,3);
a1=-A_controlable(4,4);
R=1.0;
q4=10000;
q3=10000;
q2=1;
q1=10;
p1=.5*(-4*a4*R+sqrt((4*a4*R)^2+16*q1*R));
p2=.5*(-4*a3*R+sqrt((4*a3*R)^2+16*q2*R));
p3=.5*(-4*a2*R+sqrt((4*a2*R)^2+16*q3*R));
p4=.5*(-3*a1*R+sqrt((3*a1*R)^2+8*q4*R));
P=diag([p1 p2 p3 p4]);
P(4,:)=[p1 p2 p3 p4]; %K=(([p1 p2 p3 p4])/(2*R))*inv(Mat_T);
K=(1/(2*R))*[0 0 0 1]*P*inv(Mat_T);
eig(Mat_A-Mat_B*K)
% break
while(i<(tiempo+1))
    estado=[p(i); p_p(i); alfa(i); omega(i)];
    u(i)=-K*estado;
    % u(i)=max(-200,u(i));
    % u(i)=min(200,u(i));
    p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))-Fricc*p_p(i));
    tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
    p_p(i+1)=p_p(i)+h*p_pp;
    p(i+1)=p(i)+h*p_p(i);
    omega(i+1)=omega(i)+h*tita_pp;
    alfa(i+1)=alfa(i)+h*omega(i);
    y_sal(i)=Mat_C*estado;
    i=i+1;
end
u(i)=-K*estado;
figure(1);hold on;
subplot(3,2,1);plot(t,p,color);grid on; title('Cart position');hold on;
subplot(3,2,2);plot(t,alfa,color);grid on;title('Angle');hold on;
subplot(3,2,3); plot(t,p_p,color);grid on;title('Cart velocity');hold on;
subplot(3,2,4);plot(t,omega,color);grid on;title('Angle velocity');hold on;
subplot(3,1,3);plot(t,u,color);grid on;title('Force action');xlabel('Time Sec.');hold on;
figure(2);hold on;
subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('Angle');ylabel('Angle velocity');hold on;
subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Cart position');ylabel('Cart velocity');hold on;    
