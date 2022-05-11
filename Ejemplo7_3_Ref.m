%Example extracted from: Sontag. Mathematical control theory 1998. Pag 104. http://www.sontaglab.org.
clc;clear all;
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5;
h=0.001;tiempo=(120/h);tita_pp=0;t=0:h:tiempo*h;
%Condiciones iniciales
alfa(1)=.15; color='g';
alfa(1)=2.97; color='b';
% alfa(1)=3.14; color='r';
ref=100;
omega(1)=0; p_p(1)=0; u(1)=0; p(1)=0; i=1;indice=0;
%Versión linealizada en el equilibrio inestable. Sontag Pp 104.
% estado=[p(i); p_p(i); alfa(i); omega(i)]
Mat_A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(long*M) g*(m+M)/(long*M) 0];
Mat_B=[0; 1/M; 0; -1/(long*M)];
Mat_C=[1 0 0 0]; %La salida monovariable es posición
% Construcción del sistema ampliado
Mat_Aa=[Mat_A zeros(4,1);-Mat_C 0];
Mat_Ba=[Mat_B;0];
% Mat_M=[Mat_B Mat_A*Mat_B Mat_A^2*Mat_B Mat_A^3*Mat_B];%Matriz Controlabilidad
Mat_M=[Mat_Ba Mat_Aa*Mat_Ba Mat_Aa^2*Mat_Ba Mat_Aa^3*Mat_Ba Mat_Aa^4*Mat_Ba];%Matriz Controlabilidad
%Cálculo del controlador por asignación de polos
%autovalores ecuacion caracteriatica de A lazo abierto
auto_val=eig(Mat_Aa);
c_ai=conv(conv(conv(conv([1 -auto_val(1)],[1 -auto_val(2)]),[1 -auto_val(3)]),[1 -auto_val(4)]),[1 -auto_val(5)]); %= poly(eig(Mat_A))
% c_ai=poly(auto_val);
Mat_W=[c_ai(5) c_ai(4) c_ai(3) c_ai(2) 1;
    c_ai(4) c_ai(3) c_ai(2) 1 0;
    c_ai(3) c_ai(2) 1 0 0;
    c_ai(2) 1 0 0 0;
    1 0 0 0 0];
Mat_T=Mat_M*Mat_W;
A_controlable=inv(Mat_T)*Mat_Aa*Mat_T; %Verificación de que T esté bien
% 1.0000e+00 1.0724e+01 3.7126e+01 4.8193e+01 3.2590e+01 1.0330e+01 alfa_ia=[1.0000e+00 1.0724e+01 3.7126e+01 4.8193e+01 3.2590e+011.0330e+01];%Ec caract deseadapol=fliplr(alfa_ia(2:6)-c_ai(2:6));
a5=-A_controlable(5,1);%a5
a4=-A_controlable(5,2);%a4
a3=-A_controlable(5,3);
a2=-A_controlable(5,4);
a1=-A_controlable(5,5);
R=100;
q5=100000;
q4=100000;
q3=10000;
q2=1000;
q1=10;
Q=diag([q1 q2 q3 q4 q5]);
p1=.5*(-4*a5*R+sqrt((4*a5*R)^2+16*q1*R));
p2=.5*(-4*a4*R+sqrt((4*a4*R)^2+16*q2*R));
p3=.5*(-4*a3*R+sqrt((4*a3*R)^2+16*q3*R));
p4=.5*(-4*a2*R+sqrt((4*a2*R)^2+16*q4*R));
p5=.5*(-3*a1*R+sqrt((3*a1*R)^2+ 8*q5*R));
%
% p1=.5*(-4*a5*R-sqrt((4*a5*R)^2+16*q1*R));
% p2=.5*(-4*a4*R-sqrt((4*a4*R)^2+16*q2*R));
% p3=.5*(-4*a3*R-sqrt((4*a3*R)^2+16*q3*R));
% p4=.5*(-4*a2*R-sqrt((4*a2*R)^2+16*q4*R));
% p5=.5*(-3*a1*R-sqrt((3*a1*R)^2+ 8*q5*R));
Ka=(([p1 p2 p3 p4 p5])/(2*R))*inv(Mat_T);
% [Ka,S,E]=lqr(Mat_Aa,Mat_Ba,Q,R);
H=[Mat_Aa -Mat_Ba*inv(R)*Mat_Ba'; -Q -Mat_Aa'];
[V,D]=eig(H);MX1X2=[];
for(ii=1:10)
    if D(ii,ii)<0
        MX1X2=[MX1X2 V(:,ii)];
    end
end
MX1=MX1X2(1:5,:); MX2=MX1X2(6:10,:);
P=real(MX2*inv(MX1));
% Ka=inv(R)*Mat_Ba'*P;
K=Ka(1:4); KI=-Ka(5);
eig(Mat_Aa-Mat_Ba*Ka)% break
% eig(Mat_A-Mat_B*K)
% break
while(i<(tiempo+1))
    estado=[p(i); p_p(i); alfa(i); omega(i)];
    psi_p=ref-Mat_C*estado;
    psi(i+1)=psi(i)+psi_p*h;
    u(i)=-K*estado+KI*psi(i+1);
    % u(i)=max(-2000,u(i));u(i)=min(2000,u(i));
    p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))-Fricc*p_p(i));
    tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
    p_p(i+1)=p_p(i)+h*p_pp;
    p(i+1)=p(i)+h*p_p(i);
    omega(i+1)=omega(i)+h*tita_pp;
    alfa(i+1)=alfa(i)+h*omega(i);
    y_sal(i)=Mat_C*estado;
    i=i+1;
end
u(i)=-K*estado+KI*psi(i);
figure(1);hold on;
subplot(3,2,1);plot(t,p,color);grid on; title('Cart position');hold on;
subplot(3,2,2);plot(t,alfa,color);grid on;title('Angle');hold on;
subplot(3,2,3); plot(t,p_p,color);grid on;title('Cart velocity');hold on;
subplot(3,2,4);plot(t,omega,color);grid on;title('Angle velocity');hold on;
subplot(3,1,3);plot(t,u,color);grid on;title('Force action');xlabel('Time Sec.');hold on;
figure(2);hold on;
subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('Angle');ylabel('Angle velocity');hold on;
subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Cart position');ylabel('Cart velocity');hold on;