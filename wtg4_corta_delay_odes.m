function dy=wtg4_corta_delay_odes(t,y)
%All dq variables are in the same reference frame DQ (w0*t)

global flag %wglobal torqq
%This test system has a grid equivalent impedance
%Base Values in LV: 34.5 kV and 100 MVA, @ 60 Hz
%Base Values in HV: 230 kV and 100 MVA, @ 60 Hz
%Inputs: uu=[vcdw_ref;Q_ref;w_ref];
 if t<7
     vcdw_ref=1.3;
 elseif t>=7 && t<7.5
     vcdw_ref=1.3*1.15;
 elseif t>=7.5 && t<8.5
     vcdw_ref=1.3*1.0;
 else
     vcdw_ref=1.3*0.9;
 end
Qrefw=0.0;

 if t<15
     wref=0.9;
 else
     wref=0.9*0.9;
 end

Tm=-0.7;

wb=120*pi;

%Generator side
is=y(1:2);
wr=y(3);
xintdg=y(4);
xintqg=y(5);
xintw=y(6);
%Inverter (generator side)
icwDQ=y(7:8);
vfwDQ=y(9:10);
vcdw=y(11);
delta_pcc1=y(12);
xintpll_pcc1=y(13);
xintdw=y(14);
xintqw=y(15);
xintvcdw=y(16);
xintQw=y(17);
xdelayg_d=y(18:19);
xdelayg_q=y(20:21);

xdelayr_d=y(22:23);
xdelayr_q=y(24:25);
if1DQ=y(26:27);

%Wind tubine model
w0g=2*pi*9.75;
wbg=w0g;
Rs=0.0038661;
Ld=0.453814;
Lq=0.768235;
phiPM=0.895975;
J=2*1.685;%2.3;
B=0.5;%0.2;

%Control gains generator
taug=10e-3;
kpdig=Ld/(wbg*taug);
kpqig=Lq/(wbg*taug);
kiig=Rs/taug;

kpw=5.433;%7.943;
kiw=7.16;%19.86;

%Delay switching
mf1=63;
mf2=63;
Ts1=1/(2*mf1*60);
Ts2=1/(2*mf2*60);

Td1=1.5*Ts1;
Td2=1.5*Ts2;

Adelay1=[-6/Td1 -12/(Td1^2)
                1 0];
Bdelay1=[1;0];            
Cdelay1=[-12/Td1 0];
Ddelay1=1;

Adelay2=[-6/Td2 -12/(Td2^2)
                1 0];
Bdelay2=[1;0];            
Cdelay2=[-12/Td2 0];
Ddelay2=1;

%Transformations
Tw=[cos(delta_pcc1) sin(delta_pcc1)
        -sin(delta_pcc1) cos(delta_pcc1)];

icw=Tw*icwDQ;
if1=Tw*if1DQ;

idrefg=0;
iqrefg=kpw*(wref-wr)+xintw;

udg=kpdig*(idrefg-is(1))+xintdg;
uqg=kpqig*(iqrefg-is(2))+xintqg;

ssdg=1/vcdw*(udg-wr*Lq*is(2));
ssqg=1/vcdw*(uqg+wr*Ld*is(1)+wr*phiPM);
sdg=Cdelay1*xdelayg_d+Ddelay1*ssdg;
sqg=Cdelay1*xdelayg_q+Ddelay1*ssqg;

mag=sqrt(sdg^2+sqg^2);
thettg=atan2(sqg,sdg);
if mag>=0.99
    mag=0.99;
end
sdqg=mag*[cos(thettg);sin(thettg)];

Ldq=diag([Ld,Lq]);
phis=Ldq*is+phiPM*[1;0];
Te=phis(1)*is(2)-phis(2)*is(1);

Q=[0 -1
    1 0];

Xt=Ldq;
    
A=vcdw*sdqg-Rs*is-wr*Q*phis;
dyg=[(Xt/wbg)\A
       1/J*(Te-Tm-B*wr)
       kiig*(idrefg-is(1))
       kiig*(iqrefg-is(2))
       kiw*(wref-wr)];

%Parameters (Grid side of wind turbine)
Ccdw=3.77;%1.15;
Rcw=0.0021;
Lcw=0.08;
Cfw=0.149;
Rfw=0.005;
Lfw=0.15;
w0=120*pi;

%Transformer (en p.u.)
Rt1=0.004;
Lt1=0.06;


%Constates del PLL

kppll=2.2;%Jala:5.2
kipll=750;%Jala=4600

%Impedancia de la red
Ssc=700e6;
Vb=230e3;
%wb=2*pi*60;

x_r=20;
Z=Vb^2/Ssc;
Zb=Vb^2/100e6;
Rg=(Z/sqrt(x_r^2+1))/Zb;
Lg=(Rg*x_r);

%Impedancia de la linea
%Arreglo vertical: [0+1j*28 0+1j*35 0+1j*42]
%Conductores: rcond=0.0140716, rho=0.07284/1e3;
long1=60;%60; %km
RTL=0.0751040047*long1/Zb;
LTL=0.1336645042e-2*w0*long1/Zb;
CTL=0.870412776e-8*w0*long1*Zb;

Req1=Rfw+RTL+Rt1+Rg;
Leq1=Lfw+LTL+Lt1+Lg;

Rt2=Rt1+Rg+RTL;
Lt2=Lt1+Lg+LTL;

%Punto de conexion del WTG y el lado rectificador del HVDC 
theta_plli=w0*t;
Ti_park=2/3*[cos(theta_plli) cos(theta_plli-2*pi/3) cos(theta_plli+2*pi/3)
    -sin(theta_plli) -sin(theta_plli-2*pi/3) -sin(theta_plli+2*pi/3)];
vin=1.02*[cos(w0*t+0.1)
                cos(w0*t+0.1-2*pi/3)
                cos(w0*t+0.1+2*pi/3)];
v_in=vin;
%v_in=abs(vin).*cos(w0*t+angle(vin));
vgridDQ=Ti_park*v_in;
difw_dt=w0/Leq1*(vfwDQ-Req1*if1DQ-1*Leq1*Q*if1DQ-vgridDQ);
vpccDQ=Rt2*if1DQ+Lt2/w0*difw_dt+Lt2*Q*if1DQ+vgridDQ;%(Rt1+Rfr-(Lt1+Lfr)*Req1/Leq1)*if1DQ+(Lt1+Lfr)/Leq1*vfwDQ+(1-(Lt1+Lfr)/Leq1)*vfrDQ;
vpcc1=Tw*vpccDQ;
wpcc1=kppll*vpcc1(2)+xintpll_pcc1;

%Constants
tau_in2=0.8e-3;

kpiw=(Lcw+Lfw)/(wb*tau_in2);
kiiw=(Rcw+Rfw)/(tau_in2);

kpvw=-1.258;
kivw=- 8.25;

kpQw=-0.7906;
kiQw=-17.5;

Qmw=-vpcc1(1)*if1(2)+vpcc1(2)*if1(1);

idrefw=kpvw*(vcdw_ref-vcdw)+xintvcdw;
iqrefw=kpQw*(Qrefw-Qmw)+xintQw;

udw=kpiw*(idrefw-icw(1))+xintdw;
uqw=kpiw*(iqrefw-icw(2))+xintqw;

ssdw=1/vcdw*(udw-wpcc1*(Lcw+Lfw)*icw(2)+vpcc1(1));%
ssqw=1/vcdw*(uqw+wpcc1*(Lcw+Lfw)*icw(1)+vpcc1(2));

sdw=Cdelay2*xdelayr_d+Ddelay2*ssdw;
sqw=Cdelay2*xdelayr_q+Ddelay2*ssqw;
maw=sqrt(sdw^2+sqw^2);
phiw=atan2(sqw,sdw);
if maw>=0.99
    maw=0.99;
end
sdqw=maw*[cos(phiw);sin(phiw)];
sdqwDQ=Tw\sdqw;
dygw=[wb/Lcw*(vcdw*sdqwDQ-Rcw*icwDQ-1*Lcw*Q*icwDQ-vfwDQ)
    wb/Cfw*(icwDQ-if1DQ-1*Cfw*Q*vfwDQ)
    wb/Ccdw*(-sdqg.'*is-sdqwDQ.'*icwDQ)
    wb*(wpcc1-w0/wb)
    kipll*vpcc1(2)
    kiiw*(idrefw-icw(1))
    kiiw*(iqrefw-icw(2))
    kivw*(vcdw_ref-vcdw)
    kiQw*(Qrefw-Qmw)
    Adelay1*xdelayg_d+Bdelay1*ssdg;
    Adelay1*xdelayg_q+Bdelay1*ssqg;
    Adelay2*xdelayr_d+Bdelay2*ssdw;
    Adelay2*xdelayr_q+Bdelay2*ssqw;
    difw_dt];%wb/Leq1*(vfwDQ-Req1*if1DQ-1*Leq1*Q*if1DQ-vfrDQ)];

dy=[dyg;dygw];    
if flag==1
    dy=[vpccDQ
        sdqg
        sdqw];
end
end