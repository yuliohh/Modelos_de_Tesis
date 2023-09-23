function dy=hvdc_inv_wtg4(t,y)
%t=0;
global flag
%This test system has a grid equivalent impedance
%Global sys: 150 kV and 300 MVA, @ 60 Hz

%Grid and DC-link parameters
 if t<5
     vcdw_ref=1.3;
     
 else
     vcdw_ref=1.3*1.1;
     
 end
 
 if t<7
     Qrefw=0;
 else
     Qrefw=0.15;
 end
wb=120*pi;
%State variables
icw=y(1:2);
vfw=y(3:4);
if1=y(5:6);
vcdw=y(7);
delta_pcc1=y(8);
xintpll_pcc1=y(9);
xintdw=y(10);
xintqw=y(11);
xintvcdw=y(12);
xintQw=y(13);
%Rectifier
vfr=y(14:15);
icr=y(16:17);
xintdr=y(18);
xintqr=y(19);
xintVd=y(20);
xintVq=y(21);
%Inverter
igg=y(22:23);
vfi=y(24:25);
ici=y(26:27);
deltai=y(28);
xintplli=y(29);
xintdi=y(30);
xintqi=y(31);
xintVi=y(32);

i_cd=y(33);
vcd1=y(34);
vcd2=y(35);

is=y(36:37);
wr=y(38);
xintdg=y(39);
xintqg=y(40);
xintw=y(41);

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
taug=10e-3;%100e-3;
kpdig=Ld/(wbg*taug);
kpqig=Lq/(wbg*taug);
kiig=Rs/taug;

kpw=5.433;%7.943;
kiw=7.16;%19.86;

if t<7
    wref=0.9;
    Tm=-0.8;
elseif t>=7 && t<14
    wref=0.9;
    Tm=-0.8;
elseif t>=14 && t<21
    wref=0.9;
    Tm=-0.8;
elseif t>=21 && t<28
    wref=0.9;
    Tm=-0.8;
elseif t>=21 && t<28
    wref=0.9;
    Tm=-0.8;
else 
    wref=0.9;
    Tm=-0.8;
end

idrefg=0;
iqrefg=kpw*(wref-wr)+xintw;

udg=kpdig*(idrefg-is(1))+xintdg;
uqg=kpqig*(iqrefg-is(2))+xintqg;

sdg=1/vcdw*(udg-wr*Lq*is(2));
sqg=1/vcdw*(uqg+wr*Ld*is(1)+wr*phiPM);

mag=sqrt(sdg^2+sqg^2);
thettg=atan2(sqg,sdg);
if mag>=0.99
    mag=0.99;
end
sdqg=mag*[cos(thettg);sin(thettg)];

L=diag([Ld,Lq]);
phis=L*is+phiPM*[1;0];
Te=phis(1)*is(2)-phis(2)*is(1);

Q=[0 -1
    1 0];

A=(L/wbg)\(vcdw*sdqg-Rs*is-wr*Q*phis);
dyg=[A
       1/J*(Te-Tm-B*wr)
       kiig*(idrefg-is(1))
       kiig*(iqrefg-is(2))
       kiw*(wref-wr)];

%Parameters (Grid side of wind turbine)
Ccdw=3.77;%1.15;
Rcw=0.0021;
Lcw=0.08;
Cfw=0.031;%0.149;
Rfw=0.005;
Lfw=0.15;
w0=120*pi;

%Transformer 
Rt1=0.004*2;
Lt1=0.06*2;


%Parametros (HVDC side)
Rcr=0.0032;
Lcr=0.105;
Cfr=0.1;
Rfr=0.005*2;%4.304;
Lfr=0.125*2;%4.304;

Ccd=4;
Rcd=0.005;
Lcd=0.19;%0.24;

Rci=0.0025;
Lci=0.095;
Cfi=0.1;
Rfi=0.008;
Lfi=0.12;

%Constates del PLL
kppll=0.4;%1.8736e-3*120*pi;
kipll=150;%0.3363*120*pi;

Req1=Rfw+Rt1+Rfr;
Leq1=Lfw+Lt1+Lfr;
T12=[cos(-delta_pcc1) -sin(-delta_pcc1)
          sin(-delta_pcc1) cos(-delta_pcc1)];
T21=[cos(delta_pcc1) -sin(delta_pcc1)
          sin(delta_pcc1) cos(delta_pcc1)];      
vpcc1=(Rt1+Rfr-(Lt1+Lfr)*Req1/Leq1)*if1+(Lt1+Lfr)/Leq1*vfw+(1-(Lt1+Lfr)/Leq1)*T12*vfr;
wpcc1=kppll*vpcc1(2)+xintpll_pcc1;

%Constants
tau_in2=1e-3;

kpiw=(Lcw+Lfw)/(wb*tau_in2);
kiiw=(Rcw+Rfw)/(tau_in2);

kpvw=-1.258;
kivw=- 8.25;

kpQw=-0.7906;
kiQw=-47.5;

Qmw=-vpcc1(1)*if1(2)+vpcc1(2)*if1(1);

idrefw=kpvw*(vcdw_ref-vcdw)+xintvcdw;
iqrefw=kpQw*(Qrefw-Qmw)+xintQw;

udw=kpiw*(idrefw-icw(1))+xintdw;
uqw=kpiw*(iqrefw-icw(2))+xintqw;

sdw=1/vcdw*(udw-wpcc1*(Lcw+Lfw)*icw(2)+vpcc1(1));%
sqw=1/vcdw*(uqw+wpcc1*(Lcw+Lfw)*icw(1)+vpcc1(2));
maw=sqrt(sdw^2+sqw^2);
phiw=atan2(sqw,sdw);
if maw>=0.99
    maw=0.99;
end
sdqw=maw*[cos(phiw);sin(phiw)];

dygw=[wb/Lcw*(vcdw*sdqw-Rcw*icw-wpcc1*Lcw*Q*icw-vfw)
    wb/Cfw*(icw-if1-wpcc1*Cfw*Q*vfw)
    wb/Leq1*(vfw-Req1*if1-wpcc1*Leq1*Q*if1-T12*vfr)%Grid Equation 1 Positions: (13:14)
    wb/Ccdw*(-sdqg.'*is-sdqw.'*icw)
    wb*(wpcc1-w0/wb)
    kipll*vpcc1(2)
    kiiw*(idrefw-icw(1))
    kiiw*(iqrefw-icw(2))
    kivw*(vcdw_ref-vcdw)
    kiQw*(Qrefw-Qmw)];

%HVDC
Rgg=0.002125;
Lgg=0.12412;
Req2=Rfi+Rgg;
Leq2=Lfi+Lgg;
theta_plli=w0*t+deltai;
vg2=1.02*[cos(w0*t)
    cos(w0*t-2*pi/3)
    cos(w0*t+2*pi/3)];

Ti=2/3*[cos(theta_plli) cos(theta_plli-2*pi/3) cos(theta_plli+2*pi/3)
    -sin(theta_plli) -sin(theta_plli-2*pi/3) -sin(theta_plli+2*pi/3)];

vgdqi=Ti*vg2;
vpcc2=(Rgg-Lgg*Req2/Leq2)*igg+Lgg/Leq2*vfi+(1-Lgg/Leq2)*vgdqi;
kpplli=0.4;%1.8736e-3*120*pi;
kiplli=150;%0.3363*120*pi;
wpcc2=kpplli*(vpcc2(2))+xintplli;

vpccr=vfr;
wpccr=1;

%Ganancias lado rectificador
tau_in=1e-3;
kpir=Lcr/(wb*tau_in);
kiir=Rcr/tau_in;
kpvd=0.45;
kivd=12.5;

%Ganancias lado inversor
tau_in2=1e-3;
kpii=(Lci+Lfi)/(wb*tau_in2);
kiii=(Rci+Rfi)/tau_in2;

kpvi=-2.7;
kivi=-2.7*25;

if t<5
    Vd_ref=1;
else
    Vd_ref=1.0;
end
Vq_ref=0;
iff1=T21*if1;
uvd=kpvd*(Vd_ref-vpccr(1))+xintVd;
uvq=kpvd*(Vq_ref-vpccr(2))+xintVq;
idrefr=-uvd+iff1(1)+wpccr*Cfr*vpccr(2);
iqrefr=-uvq+iff1(2)-wpccr*Cfr*vpccr(1);
edr=idrefr-icr(1);
eqr=iqrefr-icr(2);
udr=kpir*edr+xintdr;
uqr=kpir*eqr+xintqr;

sdr=1/vcd1*(-udr+wpccr*Lcr*icr(2)+vpccr(1));
sqr=1/vcd1*(-uqr-wpccr*Lcr*icr(1)+vpccr(2));%

mar=sqrt(sdr^2+sqr^2);
thettr=atan2(sqr,sdr);
if mar>0.99
  mar=0.99;
end
sdqr=mar*[cos(thettr);sin(thettr)];

Qq=[0 -1
        1 0];
dy1=[wb/Cfr*(T21*if1-icr-wpccr*Cfr*Qq*vfr)
       wb/Lcr*(vfr-Rcr*icr-wpccr*Lcr*Qq*icr-vcd1*sdqr)
       kiir*edr
       kiir*eqr
       kivd*(Vd_ref-vpccr(1))
       kivd*(Vq_ref-vpccr(2))];
   
%Estacion inversora  
if t<4
     Vrefi=1.4;
elseif t>=4 && t<10
   Vrefi=1.4*0.9;
else
     Vrefi=1.4*0.9;
end
eV=Vrefi-vcd2;
idrefi=kpvi*eV+xintVi;
iqrefi=0;
edi=idrefi-ici(1);
eqi=iqrefi-ici(2);
udi=kpii*edi+xintdi;
uqi=kpii*eqi+xintqi;


sdi=1/vcd2*(udi-wpcc2*(Lci+Lfi)*ici(2)+vpcc2(1));
sqi=1/vcd2*(uqi+wpcc2*(Lci+Lfi)*ici(1)+vpcc2(2));%

mai=sqrt(sdi^2+sqi^2);
thetti=atan2(sqi,sdi);
if mai>0.99
  mai=0.99;
end
sdqi=mai*[cos(thetti);sin(thetti)];

dy2=[wb/Leq2*(vfi-Req2*igg-wpcc2*Leq2*Qq*igg-vgdqi)%Grid Equation 2: Positions(1:2)
       wb/Cfi*(ici-igg-wpcc2*Cfi*Qq*vfi)
       wb/Lci*(vcd2*sdqi-Rci*ici-wpcc2*Lci*Qq*ici-vfi)
       wb*(wpcc2-w0/wb)
       kiplli*(vpcc2(2))
       kiii*edi
       kiii*eqi
       kivi*eV]; 

dy_hvdc=[dy1
        dy2
        wb/Lcd*(vcd1-Rcd*i_cd-vcd2)
        wb/Ccd*(sdqr.'*icr-i_cd)
        wb/Ccd*(i_cd-sdqi.'*ici)
        dyg];
        %kiv*(Vref-vpcc1(1))];

dy=[dygw;dy_hvdc];    
if flag==1
    dy=[vpcc1
            vpcc2];
end
end