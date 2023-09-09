function dy=hvdc_pu1_abcc2(t,y)
%u: Vector of references: [P,Q,Vcd]
u=[1;0;1.4];
%flag=0;
%t=0;
%Base parameters: 500 MVA, 230 kV @ 60 Hz

%Parameters
Ccd=4;
Rcd=0.005;
Lcd=0.32;

Rcr=0.012;
Lcr=0.22;
Cfr=0.15;
Rfr=0.01;
Lfr=0.16;
Rgr=0.012;
Lgr=0.22;
Rtr=Rfr+Rgr;
Ltr=Lfr+Lgr;

Rci=0.015;
Lci=0.22;
Cfi=0.17;
Rfi=0.008;
Lfi=0.12;
Rgi=0.0125;
Lgi=0.2412;
Rti=Rfi+Rgi;
Lti=Lfi+Lgi;

w0=120*pi;
wb=w0;
thetr=15*pi/180;
vgr=1.02*[cos(w0*t+thetr)
    cos(w0*t+thetr-2*pi/3)
    cos(w0*t+thetr+2*pi/3)];

theti=10*pi/180;
vgi=0.96*[cos(w0*t+theti)
    cos(w0*t+theti-2*pi/3)
    cos(w0*t+theti+2*pi/3)];

%State-variables
igr=y(1:3);
vfr=y(4:6);
icr=y(7:9);
deltar=y(10);
xintpllr=y(11);
xintdr=y(12);
xintqr=y(13);
xintPr=y(14);
xintQr=y(15);

igi=y(16:18);
vfi=y(19:21);
ici=y(22:24);
deltai=y(25);
xintplli=y(26);
xintdi=y(27);
xintqi=y(28);
xintvi=y(29);

vcdr=y(30);
icd=y(31);
vcdi=y(32);

theta_pllr=w0*t+deltar;
kppllr=1.8736e-3*120*pi;
kipllr=0.3363*120*pi;

Tr=(2/3)*[cos(theta_pllr) cos(theta_pllr-2*pi/3) cos(theta_pllr+2*pi/3)
    -sin(theta_pllr) -sin(theta_pllr-2*pi/3) -sin(theta_pllr+2*pi/3)];

theta_plli=w0*t+deltai;
kpplli=1.8736e-3*120*pi;
kiplli=0.3363*120*pi;

Ti=(2/3)*[cos(theta_plli) cos(theta_plli-2*pi/3) cos(theta_plli+2*pi/3)
    -sin(theta_plli) -sin(theta_plli-2*pi/3) -sin(theta_plli+2*pi/3)];

%Voltage in PCC (rectifier side)
vpccr=Lgr/Ltr*((Rfr-Rgr*Lfr/Lgr)*igr+Lfr/Lgr*vgr+vfr);
vpccr0=Tr*vpccr;
icr0=Tr*icr;
igr0=Tr*igr;
wr=kppllr*vpccr0(2)+xintpllr;

%Voltage in PCC (inverter side)
vpcci=Lfi/Lti*((Rgi-Rfi*Lgi/Lfi)*igi+Lgi/Lfi*vfi+vgi);
vpcci0=Ti*vpcci;
ici0=Ti*ici;
wi=kpplli*vpcci0(2)+xintplli;

%Rectifier Side
taur=1e-3;
kpir=(Lcr+Lfr)/(taur*wb);
kiir=(Rcr+Rfr)/taur;

kppr=0.435;
kipr=130;

kpqr=-0.435;
kiqr=-130;

%Reference chanages
if t<0.5
    Prefr=u(1);
    Qrefr=u(2);
else
   Prefr=u(1)*1.15;
   Qrefr=0.1;
end

Pmeasr=vpccr0(1)*igr0(1)+vpccr0(2)*igr0(2);
errPr=Prefr-Pmeasr;
idrefr=kppr*errPr+xintPr;

Qmeasr=-vpccr0(1)*igr0(2)+vpccr0(2)*igr0(1);
errQr=Qrefr-Qmeasr;
iqrefr=kpqr*errQr+xintQr;

errdr=idrefr-icr0(1);
errqr=iqrefr-icr0(2);
udr=kpir*errdr+xintdr;
uqr=kpir*errqr+xintqr;

sdr=1/vcdr*(-udr+wr*(Lcr+Lfr)*icr0(2)+vpccr0(1));
sqr=1/vcdr*(-uqr-wr*(Lcr+Lfr)*icr0(1)+vpccr0(2));

mar=sqrt(sdr^2+sqr^2);
thettr=atan2(sqr,sdr);

 if mar>0.99
     mar=0.99;
 end
sdqr=mar*[cos(thettr);sin(thettr)];
mabcr=(3/2*Tr).'*sdqr;
% %CONMUTADO
mf1=27;
tri1=sawtooth(mf1*w0*t+1.5*pi,0.5);

sar=mabcr(1)>tri1;%1/2*(1+tanh(mf1*(mabcr(1)-tri1)));
sbr=mabcr(2)>tri1;%1/2*(1+tanh(mf1*(mabcr(2)-tri1)));%
scr=mabcr(3)>tri1;%1/2*(1+tanh(mf1*(mabcr(3)-tri1)));
sabcr=[sar;sbr;scr];

ssumr=1/3*sum(sabcr);
g1r=sabcr-ssumr;
%Set of ODEs in rectifier side
dyr=[wb/Ltr*(vgr-Rtr*igr-vfr)
       wb/Cfr*(igr-icr)
       wb/Lcr*(vfr-Rcr*icr-2*vcdr*g1r)
       wb*(wr-1)
       kipllr*vpccr0(2)
       kiir*errdr
       kiir*errqr
       kipr*errPr
       kiqr*errQr];
 %%%----------------------------------------------------------------  
%Inverter Side  
taui=1e-3;
kpii=(Lci+Lfi)/(taui*wb);
kiii=(Rci+Rfi)/taui;

kpvi=-2.7;
kivi=-2.7*25;

%Cahnge of reference
 if t<0.2
    vcd_refi=u(3);
 else
    vcd_refi=1.55;
 end

errvi=vcd_refi-vcdi;
idrefi=kpvi*errvi+xintvi;

iqrefi=0;

errdi=idrefi-ici0(1);
errqi=iqrefi-ici0(2);
udi=kpii*errdi+xintdi;
uqi=kpii*errqi+xintqi;

sdi=1/vcdi*(udi-wi*(Lci+Lfi)*ici0(2)+vpcci0(1));
sqi=1/vcdi*(uqi+wi*(Lci+Lfi)*ici0(1)+vpcci0(2));

mai=sqrt(sdi^2+sqi^2);
thetti=atan2(sqi,sdi);

 if mai>0.99
     mai=0.99;
 end
sdqi=mai*[cos(thetti);sin(thetti)];
mabci=(3/2*Ti).'*sdqi;
% %CONMUTADO
mf2=27;
tri2=sawtooth(mf2*w0*t+1.5*pi,0.5);
sai=mabci(1)>tri2;%1/2*(1+tanh(mf2*(mabci(1)-tri2)));%
sbi=mabci(2)>tri2;%1/2*(1+tanh(mf2*(mabci(2)-tri2)));%
sci=mabci(3)>tri2;%1/2*(1+tanh(mf2*(mabci(3)-tri2)));%
sabci=[sai;sbi;sci];

ssumi=1/3*sum(sabci);
g2i=sabci-ssumi;
%Set of ODEs in inverter station
dyi=[wb/Lti*(-vgi-Rti*igi+vfi)
       wb/Cfi*(-igi+ici)
       wb/Lci*(-vfi-Rci*ici+2*vcdi*g2i)
       wb*(wi-1)
       kiplli*vpcci0(2)
       kiii*errdi
       kiii*errqi
       kivi*errvi];   
%System dx/dt=f(x,t)   
dy=[dyr
        dyi
        wb/Ccd*((4/3*sabcr).'*icr-icd)
        wb/Lcd*(vcdr-Rcd*icd-vcdi)
        wb/Ccd*(icd-(4/3*sabci).'*ici)];
end