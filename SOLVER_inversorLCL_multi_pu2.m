%SOLVER inversorLCL_multi_pu numerico
%clc
clear all
global Ri Li Rd Cf R L Rf Lf KP KI VCD IDREF IQREF MF Lg Rg 
global Ccap kp_pll ki_pll Vp f num h
global Vb Sb Ib Zb Anb wb Xb Lb Cb
%Parametros de las n unidades vsc
num=2;
Rg=3.85278e-3;%0.65813e-3*2.20;
Lg=0.102198e-3;%9.8969e-6*0.5;

%vsc1
Ri1=0.64533e-3;%0.1;
Li1=0.435e-3;%1.2e-3;

Rd1=6.3e-3;%2;
Cf1=5.3e-3;%10e-6;

R1=0.32266e-3;%0.1;
L1=1.33e-6;%1e-3;

Rf1=0.2*Rg;
Lf1=0.2*Lg;

%vsc2
Ri2=0.968e-3;%0.1;
Li2=0.5135e-3;%1.2e-3;

Rd2=0.2484e-3;%2;
Cf2=8.7e-3;%10e-6;

R2=0.70986e-3;%0.1;
L2=0.6012e-6;%1e-3;

Rf2=0.2*Rg;
Lf2=0.2*Lg;

%red electrica
Ccap=120e-6;
Rt1=Ri1+R1;
Lt1=Li1+L1;

Rt2=Ri2+R2;
Lt2=Li2+L2;

kp_pll=4;
ki_pll=4e3;

mf1=27;
mf2=21;
f=60;
Vp=480*sqrt(2/3);
Vcd1=1800;
Vcd2=1800;

Idref1=2*0.9*1e6/(3*Vp);
Iqref1=0;

Idref2=2*0.85*1.5e6/(3*Vp);
Iqref2=0;

h=45;
tau=0.01e-3;
kp1=4*Lt1/tau;
ki1=4*Rt1/tau;

kp2=4*Lt2/tau;
ki2=4*Rt2/tau;
PerD=10/100;

%Valores base
Vb=Vp/sqrt(2);
Sb=1.5e6/3;
Ib=Sb/Vb;
Zb=Vb^2/Sb;
Anb=1;
wb=2*pi*f;
Lb=Zb/wb;
Cb=1/(Zb*wb);
Xb=10e5;

%Vectores con las cantidades de cada vsc en pu's
Ri=[Ri1 Ri2]/Zb;
Li=[Li1 Li2]/Lb;
Rd=[Rd1 Rd2]/Zb;
Cf=[Cf1 Cf2]/Cb;
R=[R1 R2]/Zb;
L=[L1 L2]/Lb;
Rf=[Rf1 Rf2]/Zb;
Lf=[Lf1 Lf2]/Lb;

%%Ganancias modificadas para ver inestabilidad
%KP=[kp1*238.7e-3 kp2];
%%Ganancias normalonas
KP=[kp1 kp2];
KI=[ki1 ki2];
VCD=[Vcd1 Vcd2]/Vb;
IDREF=[Idref1 Idref2]/Ib;
IQREF=[Iqref1 Iqref2]/Ib;
MF=[mf1 mf2];
N=5001;
u=[IDREF(1);IQREF(1);IDREF(2);IQREF(2)];
opt=odeset('reltol',1e-6,'abstol',1e-6);
[t,y]=ode23tb('inversorLCL_multi_pu2',[0 3/f],zeros(13*num+6,1),opt);
plot(t,y(:,1:3))
