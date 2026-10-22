% num=[45 1]
% PRBS 7:  7 6
% PRBS 9:  9 5
% PRBS 11:  11 9
% PRBS 13: 13 12 2 1
% PRBS 15:  15 14
% PRBS 20:  20 3
% PRBS 23:  23 18
% PRBS 31: 31 28
% num=1;
% den= ([0.064 .5024 .8 1]);
% Éste algoritmo funciona porque se ha ajustado a la señal de entrada y al
% tiempo de muestro a valores adecuados.
clear all; close all;clc;
Med=2500;
% ts=1/150; orden_a=3; orden_b=2;  num=1;%Es tercer orden, 3-2, aprox 2-1 para num=1;
ts=1/200; orden_a=3; orden_b=3; num=[45 1];%Es tercer orden, 3-2, aprox 2-1 para num=[45 1];
t=0:ts:20000*ts;
randn('state',0)
% Señal de entrada
StepAmplitude=1.5;
ue=StepAmplitude*sign(sin(2*pi*.2*1.0*t));
m=31;
x=ones(m,1);
N=length(ue);%Puntos de la PRBS para muestrear
el=50;%Repeticiones del ZOH
for k=1:el:N
    n_b=xor(x(18),x(23));
    y(k:k+el-1)=x(7);
    for h=m-1:-1:1
        x(h+1)=x(h);
    end
    x(1)=n_b; %Ingreso el nuevo valor
end
% x=[];
% plot(y)
% break
ue=y(1:N); %Uso de PRBS
%Planta Ec 8-23
den= [0.064 .5024 .8 1];
sys=tf(num,den);
[y_D,t_D]=lsim(sys,ue,t);
% plot(y_D);hold on;%Respuesta
% break
ys=y_D';
%Cuanto es lo que tarda en iniciar el algoritmo
off_set=orden_a+1+50;
u=zeros(Med,1);
z=zeros(Med,1);
u=ue(1+off_set:length(u)+off_set)';ui=u;
z=ys(1+off_set:length(z)+off_set)';zi=z;
for jj=orden_a+1:Med
    vec_a=fliplr([ u(jj-orden_b:jj-1); -z(jj-orden_a:jj-1)]');
    H(jj-orden_a,:)=vec_a;
end
[aa bb ]=size(H);
Z=(z(orden_a+1:end));
in_1=H';
in_2=in_1*H;
in_3=inv(in_2);
in_4=in_3*in_1;
c=in_4*Z; % Ec 131
abs(roots([1; c(1:orden_a)]))
zo=z;
z=zeros(Med,1);u=ui;
z=zi(1:orden_a);
u=ue(1+off_set:off_set+Med)';
for k=orden_a+1:length(u)
    zt=-flip(z(k-orden_a:k-1));
    ut=flip(u(k-orden_b:k-1));
    z(k)=c'*[zt;ut];    %Ec 125
end
dend=[1; c(1:orden_a)]';numd=[c(orden_a+1:end)]';
sys_id=tf(numd,dend,ts,'Variable','z^-1');
ue=sign(sin(2*pi*.010*t));
[y_sal,t_sal]=lsim(sys_id,ue,t,[0,0]);
[y_D,t_D]=lsim(sys,ue,t,[0,0]);
subplot(2,1,1);hold on;
plot(z,'.');
plot(zo,'+r');
plot(ys(1+off_set:off_set+Med),'y');legend('Estimada','Mediciones','Real');legend('Boxoff'),
title(['Ajuste con el orden b_s/a_s:' num2str(orden_b) '/' num2str(orden_a)]);
xlabel('Muestras')
subplot(2,1,2)
hold on
plot(t_D*ts,y_D,'.');
plot(t_sal*ts,y_sal,'k');legend('Real','Identificada'),legend('Boxoff')
title('Desempeño del modelo ajustado');xlabel('Tiempo. [Seg.]')
figure
step(StepAmplitude*sys,'r',StepAmplitude*sys_id,'k'),hold on
legend('Real','Identificada');
legend('boxoff');
sys %Sistema original
sysc = d2c(sys_id,'tustin') %sistema identificado
figure
bode(sys,logspace(-1,1));
hold on;
bode(sys_id,logspace(-1,1));
legend('Original','Identificada');
bode(sysc,logspace(-1,1));
legend('boxoff')
