function dy=Two_wtg4_delay_lumped_v4(t,y)
%All dq variables are in the same reference frame DQ (w0*t)

%PMSG: without damping windings (like simulink model)
global flag wglobal torqq
%This test system has two WTG-4 interconnected 
%by means transmission lines to a grid equivalent
%Base system: 100 MVA, 34.5 kV/230 kV, 60 Hz

%Inputs: uu=[vcdw_ref;Q_ref;w_ref];
 if t<6
      vcdw_ref1=1.3;
 elseif t>=6 && t<6.5
     vcdw_ref1=1.3;
 else
     vcdw_ref1=1.3;
 end
 
Qrefw1=0.0;

% if t<10
     wref1=wglobal;%
% elseif t>=10 && t<20
%     wref1=wglobal+0.1;%
% else
%     wref1=wglobal;%
% end
Tm1=-torqq;%

%Generator 1 side
is1=y(1:2);
wr1=y(3);
xintdg1=y(4);
xintqg1=y(5);
xintw1=y(6);
%Inverter 1 (generator side)
icwDQ1=y(7:8);
vfwDQ1=y(9:10);
vcdw1=y(11);
delta_pcc1=y(12);
xintpll_pcc1=y(13);
xintdw1=y(14);
xintqw1=y(15);
xintvcdw1=y(16);
xintQw1=y(17);
if1DQ1=y(18:19);

v_pi1=y(20:21);
i_tl1=y(22:23);
v_picom=y(24:25);


igridDQ=y(26:27);

w0=120*pi;
%Transmission Line 1 
Zb=230^2/100;
long1=130;%60; %km
R_tl1=0.0751040047*long1/Zb;
X_tl1=0.1336645042e-2*w0*long1/Zb;
Y_tl1=0.870412776e-8*w0*long1*Zb;

%Grid equivalent
Ssc=1500e6;
Vb=230e3;
Sb=100e6;
wb=2*pi*60;
Zb=Vb^2/Sb;
Lb=Zb/wb;

x_r=20;
Z=Vb^2/Ssc;
Rr=Z/sqrt(x_r^2+1);
Lr=Rr*x_r/wb;

Rgrid= Rr/Zb;
Lgrid=Lr/Lb;

%%
%Wind tubine model 1
w0g1=2*pi*9.75;%(Electrical and poles=48)
wbg1=w0g1;
Rs1=0.0038661;
Ld1=0.453814;
Lq1=0.768235;
phiPM1=0.895975;
J1=2*1.685;%2.3;
B1=0.5;%0.2;

%Control gains generator 1
taug=10e-3;
kpdig1=Ld1/(wbg1*taug);
kpqig1=Lq1/(wbg1*taug);
kiig1=Rs1/taug;

kpw1=5.433;
kiw1=7.16;

%Transformations
Tw1=[cos(delta_pcc1) sin(delta_pcc1)
        -sin(delta_pcc1) cos(delta_pcc1)];

icw1=Tw1*icwDQ1;
if1=Tw1*if1DQ1;

idrefg1=0;
iqrefg1=kpw1*(wref1-wr1)+xintw1;

udg1=kpdig1*(idrefg1-is1(1))+xintdg1;
uqg1=kpqig1*(iqrefg1-is1(2))+xintqg1;

sdg1=1/vcdw1*(udg1-wr1*Lq1*is1(2));
sqg1=1/vcdw1*(uqg1+wr1*Ld1*is1(1)+wr1*phiPM1);

mag1=sqrt(sdg1^2+sqg1^2);
thettg1=atan2(sqg1,sdg1);
if mag1>=0.99
    mag1=0.99;
end
sdqg1=mag1*[cos(thettg1);sin(thettg1)];

Ldq1=diag([Ld1,Lq1]);
phis1=Ldq1*is1+phiPM1*[1;0];
Te1=phis1(1)*is1(2)-phis1(2)*is1(1);

Q=[0 -1
    1 0];

Xt1=Ldq1;
    
A1=vcdw1*sdqg1-Rs1*is1-wr1*Q*phis1;
dyg1=[(Xt1/wbg1)\A1
       1/J1*(Te1-Tm1-B1*wr1)
       kiig1*(idrefg1-is1(1))
       kiig1*(iqrefg1-is1(2))
       kiw1*(wref1-wr1)];

%Parameters (Grid side of wind turbine 1)
Ccdw1=8;%1.15;
Rcw1=0.008;
Lcw1=0.18;
Cfw1=0.15;
Rfw1=0.006;
Lfw1=0.11;


%Transformer 1
zz1=0.1;
x_r=25;
Rt1=zz1/sqrt(x_r^2+1);
Lt1=x_r*Rt1;


%Constates del PLL 1
kppll1=0.25;
kipll1=81.62;


Req1=Rfw1+Rt1;
Leq1=Lfw1+Lt1;

%Punto de conexion del WTG y el lado rectificador del HVDC 
theta_plli=w0*t;
Ti_park=2/3*[cos(theta_plli) cos(theta_plli-2*pi/3) cos(theta_plli+2*pi/3)
    -sin(theta_plli) -sin(theta_plli-2*pi/3) -sin(theta_plli+2*pi/3)];
vin=1.02*[cos(w0*t+0.1)
                cos(w0*t+0.1-2*pi/3)
                cos(w0*t+0.1+2*pi/3)];

v_in=vin;%abs(vin).*cos(w0*t+angle(vin));
vgridDQ=Ti_park*v_in;
difw_dt1=w0/Leq1*(vfwDQ1-Req1*if1DQ1-1*Leq1*Q*if1DQ1-v_pi1);
vpccDQ1=Rt1*if1DQ1+Lt1/w0*difw_dt1+Lt1*Q*if1DQ1+v_pi1;%(Rt1+Rfr-(Lt1+Lfr)*Req1/Leq1)*if1DQ+(Lt1+Lfr)/Leq1*vfwDQ+(1-(Lt1+Lfr)/Leq1)*vfrDQ;
vpcc1=Tw1*vpccDQ1;
wpcc1=kppll1*vpcc1(2)+xintpll_pcc1;

%Constants
tau_in21=0.8e-3;

kpiw1=(Lcw1+Lfw1)/(wb*tau_in21);
kiiw1=(Rcw1+Rfw1)/(tau_in21);

kpvw1=-2.182;
kivw1=- 28.37;

kpQw1=-0.1;
kiQw1=-2.5;

Qmw1=-vpcc1(1)*if1(2)+vpcc1(2)*if1(1);

idrefw1=kpvw1*(vcdw_ref1-vcdw1)+xintvcdw1;
iqrefw1=kpQw1*(Qrefw1-Qmw1)+xintQw1;

udw1=kpiw1*(idrefw1-icw1(1))+xintdw1;
uqw1=kpiw1*(iqrefw1-icw1(2))+xintqw1;

sdw1=1/vcdw1*(udw1-1*(Lcw1+Lfw1)*icw1(2)+vpcc1(1));%
sqw1=1/vcdw1*(uqw1+1*(Lcw1+Lfw1)*icw1(1)+vpcc1(2));
maw1=sqrt(sdw1^2+sqw1^2);
phiw1=atan2(sqw1,sdw1);
if maw1>=0.99
    maw1=0.99;
end
sdqw1=maw1*[cos(phiw1);sin(phiw1)];
sdqwDQ1=Tw1\sdqw1;
dygw1=[wb/Lcw1*(vcdw1*sdqwDQ1-Rcw1*icwDQ1-1*Lcw1*Q*icwDQ1-vfwDQ1)
    wb/Cfw1*(icwDQ1-if1DQ1-1*Cfw1*Q*vfwDQ1)
    wb/Ccdw1*(-sdqg1.'*is1-sdqwDQ1.'*icwDQ1)
    wb*(wpcc1-w0/wb)
    kipll1*vpcc1(2)
    kiiw1*(idrefw1-icw1(1))
    kiiw1*(iqrefw1-icw1(2))
    kivw1*(vcdw_ref1-vcdw1)
    kiQw1*(Qrefw1-Qmw1)
    difw_dt1];%wb/Leq1*(vfwDQ-Req1*if1DQ-1*Leq1*Q*if1DQ-vfrDQ)];

Y_eqcom=Y_tl1/2;
%Transmission Line 1
dy_tl1=[wb/(Y_tl1/2)*(if1DQ1-i_tl1-1*(Y_tl1/2)*Q*v_pi1)
            wb/X_tl1*(v_pi1-R_tl1*i_tl1-1*X_tl1*Q*i_tl1-v_picom)
            wb/Y_eqcom*(i_tl1-igridDQ-1*Y_eqcom*Q*v_picom)];
%%
%Grid Equivalent   
dy_grid=w0/Lgrid*(v_picom-Rgrid*igridDQ-1*Lgrid*Q*igridDQ-vgridDQ);
%%
dy=[dyg1;dygw1;dy_tl1;dy_grid];    
if flag==1
    dy=[vpccDQ1
        sdqg1
        sdqw1];
end
end