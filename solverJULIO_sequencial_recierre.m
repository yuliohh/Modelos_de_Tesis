%Solucion del problema para cierre sequencial de la linea, con los tiempos
clc
clear all
close all
%Sequential closer data
ta=3e-3;
tb=5e-3;
tc=7e-3;
rsw=2.5e-3;%Resistencia de breaker
%DEFINE TRANSMISSION LINE LENGTH
tic
long = 100e3; %length

%simulation parameters
Tobs = 1/60; %simulation time
N    = 2^12; %Number of samples
%Tobs = Tobs/0.95;
K = floor(0.95*N);
m     = [1:2:2*N];
delt = Tobs/N;
delw  = pi/Tobs;
n    = [0:N-1];
t    = n*delt;
c    = 2*delw;
s    = c + j*m*delw;

%Study cases


%Select the configuration's line 
rho_earth=100;
%Select the configuration's line (solo se correra el caso n=5 (LT=3) y n=7 (LT=5))
LT=5;

if LT == 1
    %HSIL Coordinates
    x = [-1.9719 -3.2902 -1.9719 -3.2902 ...
        0.0867 -0.0867  0.0867 -0.0867 ...
        1.9719  3.2902  1.9719  3.2902];
    y = [32.6975 32.9918 32.3025 32.0082 ...
        37.7857 37.7857 37.0870 37.0870 ...
        32.6975 32.9918 32.3025 32.0082];
    Pos=x+j*y;    
    Ncond=12;
    Ncpp=4;
    %HSIL conductor type and DC resistance
    rcond = 1.14*0.0254/2;
    Rcd = 0.05962/1000; %catbird
elseif LT==2
    %coordinates of conventional line
    x = [-10.7715 -11.2285 -11.2285 -10.7715...
        0.2285  -0.2285  -0.2285   0.2285...
        11.2285  10.7715  10.7715  11.2285];
    y = [28.4285 28.4285 27.9715 27.9715...
        28.4285 28.4285 27.9715 27.9715...
        28.4285 28.4285 27.9715 27.9715];
    %coodinates of equivalent conductors
    x=[-11 0 11];
    y=[28.2 28.2 28.2];
    Pos=x+j*y;    
    Ncond=3;
    Ncpp=1;
    %conductor type and DC resistance of conventional line
    rcond = 1.165*0.0254/2;
    %rcond=(4*(0.64629/2)^(4-1)*rcond)^(1/4);%Formula especifica para haz de 4 conductores de dimetro 0.64629 m
    Rcd = 0.05992/1000; %rail
    D_haz=0.64629;
elseif LT==3
    %5 conductors per phase
    Pos=[-2.5017 + 33.2785038854688i,-3.63214804288976 + 33.0186093069642i,-2.86485474077037 + 31.3420421025726i,-2.14318705375854 + 31.8344902504318i,-2.02831875234918 + 32.9042878626047i...
    0.211602690825290 + 37.2079756772756i,0.134308930148830 + 36.7318691760860i,-0.211602690825290 + 37.2079756772756i,-0.134308930148830 + 36.7318691760860i,1.39915896802585e-17 + 37.1452295593006i...
    2.50171231720211 + 33.2785038854688i,3.63214804288976 + 33.0186093069642i,2.86485474077037 + 31.3420421025726i,2.14318705375854 + 31.8344902504318i,2.02831875234918 + 32.9042878626047i];
    x=real(Pos)';
    y=imag(Pos)';
    Ncond=15;
    Ncpp=5;
    rcond=1.036*0.0254/2;
    Rcd=0.0795/1e3;%stilt
elseif LT==4    
    %6 conductors per phase
    Pos=[-1.9783 + 33.0077i,-2.8784 + 33.2623i,-1.9783 + 31.9923i,-2.8784 + 31.7377i,-1.6800 + 32.5000i,-3.1200 + 32.5000i...
    0.0823 + 36.9452i,-0.0830 + 36.9452i,0.0830 + 36.3686i,-0.0830 + 36.3686i,0.2285 + 36.6569i,-0.2285 + 36.6569i...
    1.9783 + 33.0077i,2.8784 + 33.2623i,1.9783 + 31.9923i,2.8784 + 31.7377i,1.7609 + 32.5000i,3.0391 + 32.5000i];
    x=real(Pos)';
    y=imag(Pos)';
    Ncond=18;
    Ncpp=6;
    rcond=0.953*0.0254/2;
    Rcd=0.10236/1e3;%eagle
elseif LT==5
    %7 conductors per phase
    Pos=[-2.0309 + 33.1380i,-2.9686 + 33.1706i,-3.3889 + 32.6626i,-2.7138 + 31.4156i,-1.9598 + 31.5705i,-1.6590 + 32.1219i,-1.7655 + 32.8465i...
    0.2299 + 36.7319i,0.2399 + 36.3189i,0.0991 + 36.1046i,-0.2299 + 36.7319i,-0.2399 + 36.3189i,-0.0991 + 36.1046i,1.3992e-17 + 36.5390i...
    2.0309 + 33.1380i,2.9686 + 33.1706i,3.3889 + 32.6626i,2.7138 + 31.4156i,1.9598 + 31.5705i,1.6590 + 32.1219i,1.7655 + 32.8465i];
    x=real(Pos)';
    y=imag(Pos)';
    Ncond=21;
    Ncpp=7;
    rcond=0.858*0.0254/2;
    Rcd=0.1195/1e3;%hawk
elseif LT==6
    %8 conductors per phase
    Pos=[-1.3003 + 32.6162i,-1.6663 + 32.9595i,-2.2683 + 33.2121i,-3.0561 + 32.9653i,-1.3003 + 32.3838i,-1.6663 + 32.0405i,-2.2683 + 31.7879i,-3.0561 + 32.0347i...
    0.0955 + 36.0977i,0.2988 + 35.9040i,-0.0955 + 36.0977i,-0.2988 + 35.9040i,0.0955 + 35.6573i,0.2988 + 35.8510i,-0.0955 + 35.6573i,-0.2988 + 35.8510i...
    1.3003 + 32.6162i,1.6663 + 32.9595i,2.2683 + 33.2121i,3.0561 + 32.9653i,1.3003 + 32.3838i,1.6663 + 32.0405i,2.2683 + 31.7879i,3.0561 + 32.0347i];
    x=real(Pos)';
    y=imag(Pos)';
    Ncond=24;
    Ncpp=8;
    rcond=0.824*0.0254/2;
    Rcd=0.1151/1e3;%tailorbird
elseif LT==0
    Pos=[-15+j*45 0+j*45 15+j*45];
    x=real(Pos)';
    y=imag(Pos)';
    Ncond=3;
    Ncpp=1;
    rcond=0.0169;
    Rcd=0.48/1e3;%tailorbird

end

%%%%%%%%%%%%%%%%%%%%%  E  X  A  M  P  L  E  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HSIL Coordinates

%Pos=[-15+j*31, j*27, 15+j*31];
%Ncond=3;
%Ncpp=1;
%%HSIL conductor type and DC resistance
%rcond = 2.91/100;
%Rcd = 0.031/1000; 
%D_haz=0.64629;

%%%%%%%%%%%%%%%%%%%%%  E  X  A  M  P  L  E  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Source (3 phase AC)
theta = 0; w0 = 2*pi*60;
ang = theta*pi/180;
ang2 = (theta-120)*pi/180;
ang3 = (theta+120)*pi/180;

va=cos(w0*t+ang);
vb=cos(w0*t+ang2);
vc=cos(w0*t+ang3);


Vs=zeros(3,N);
Vs(1,:) = DNLT(va,c,t(2));
Vs(2,:) = DNLT(vb,c,t(2));
Vs(3,:) = DNLT(vc,c,t(2));

%Admitance source matrix
Gs=eye(Ncond)*1e6;

Na=find(t>=ta);
Na=Na(1);
maskA=ones(1,N);
maskA(1:Na)=0;

Nb=find(t>=tb);
Nb=Nb(1);
maskB=ones(1,N);
maskB(1:Nb)=0;

Nc=find(t>=tc);
Nc=Nc(1);
maskC=ones(1,N);
maskC(1:Nc)=0;

ISS=zeros(Ncond,N);
Nph=Ncond/3;
for k2=1:3
  for k=1:Nph
    ISS(k+Nph*(k2-1),:)=Vs(k2,:)*1e6;
  end
end

%Current injection vector
IT = zeros(3*Ncond,N);

IT(1:Ncond,:)=ISS;

%Load Admitance
YL=zeros(Ncond);%eye(Ncond)*3e3;
k0=[1:Ncond];
for k1=1:3
  kk=k0((k1-1)*Ncpp+1:(k1)*Ncpp);
  for kk1=1:Ncpp
    ll=kk(kk1);
    for kk2=1:Ncpp

    ll2=kk(kk2);
      if ll~=ll2
        YL(ll,ll2)=-1e6;
      else
        YL(ll,ll2)=1e6*(Ncpp-1);
      end
    end 
  end
end

% Inverse NLT for transformation of V to time domain (odd sampling)
sigma = 0.5*(1 + cos(0.5*pi*m/N)); %window function (Hanning)
cte  = 2*N*delw/pi;
exp1 = exp((c*delt+1i*pi/N).*n);
Cn   = cte*exp1;

%Loop for frequency sweep matrices A and B
for k=1:N
    [Z,Y]=TL_param_Nhaz_JULIO(s(k),Pos,[rcond Rcd Ncond rho_earth 0 1 0]);
    [M,val] = eig(Z*Y); %Modal decomposition
    cprop = sqrt(val);
    Y0 = Z\M*diag(diag(cprop))/M; %characteristic admittance matrix
    x2 = cprop.*long;
    %Admittance model of transmission line
    A(:,:,k) = Y0*M*(diag(coth(diag(x2))))/M;
    B(:,:,k) = Y0*M*(diag(csch(diag(x2))))/M;
end
O=zeros(Ncond);

%Condicion antes del cierre de los interruptores
Gsw=zeros(Ncond);
for k = 1:N
    Ybus = [Gsw+Gs -Gsw O
            -Gsw A(:,:,k)+Gsw -B(:,:,k)
            O -B(:,:,k) A(:,:,k)+YL];
    V1(:,k) = Ybus\IT(:,k); 
end

Vf=V1;
for xx = 1:3*Ncond
    fn  = ifft(Vf(xx,:).*(sigma));
    vt(xx,:)   = real(Cn.*fn);
end


%Condicion para el cierre de la fase A en t=t_a ms
ISW=zeros(Ncond,N);
IT=zeros(3*Ncond,N);

for k=1:Nph
  vf=(va-vt(Ncond+k,:)).*maskA;
  ISW(k,:)=DNLT(vf,c,t(2))/rsw;
end

IT(1:Ncond,:)=-ISW;
IT(Ncond+1:2*Ncond,:)=ISW;

Gsw=diag([1/rsw*ones(1,Nph) zeros(1,2*Nph)]);
for k = 1:N
    Ybus = [Gsw+Gs -Gsw O
            -Gsw A(:,:,k)+Gsw -B(:,:,k)
            O -B(:,:,k) A(:,:,k)+YL];
    V2(:,k) = Ybus\IT(:,k);
end

Vf=V1+V2;
for xx = 1:3*Ncond
    fn  = ifft(Vf(xx,:).*(sigma));
    vt(xx,:)   = real(Cn.*fn);
end

%Condicion para el cierre de la fase B en t=t_b ms
ISW=zeros(Ncond,N);
IT=zeros(3*Ncond,N);

for k=1:Nph
  vf=(vb-vt(Ncond+Nph+k,:)).*maskB;
  ISW(Nph+k,:)=DNLT(vf,c,t(2))/rsw;
end

IT(1:Ncond,:)=-ISW;
IT(Ncond+1:2*Ncond,:)=ISW;

Gsw=diag([1/rsw*ones(1,2*Nph) zeros(1,Nph)]);
for k = 1:N
    Ybus = [Gsw+Gs -Gsw O
            -Gsw A(:,:,k)+Gsw -B(:,:,k)
            O -B(:,:,k) A(:,:,k)+YL];
    V3(:,k) = Ybus\IT(:,k);
end

Vf=V1+V2+V3;
for xx = 1:3*Ncond
    fn  = ifft(Vf(xx,:).*(sigma));
    vt(xx,:)   = real(Cn.*fn);
end


%Condicion para el cierre de la fase C en t=t_c ms
ISW=zeros(Ncond,N);
IT=zeros(3*Ncond,N);

for k=1:Nph
  vf=(vc-vt(Ncond+2*Nph+k,:)).*maskC;
  ISW(2*Nph+k,:)=DNLT(vf,c,t(2))/rsw;
end

IT(1:Ncond,:)=-ISW;
IT(Ncond+1:2*Ncond,:)=ISW;

Gsw=eye(Ncond)*(1/rsw);
for k = 1:N
    Ybus = [Gsw+Gs -Gsw O
            -Gsw A(:,:,k)+Gsw -B(:,:,k)
            O -B(:,:,k) A(:,:,k)+YL];
    V4(:,k) = Ybus\IT(:,k);
end

Vf=V1+V2+V3+V4;
for xx = 1:3*Ncond
    fn  = ifft(Vf(xx,:).*(sigma));
    vt(xx,:)   = real(Cn.*fn);
end
toc
vt=vt';
t=t';
v0=[t vt];
v=1;
plot(v0(:,1),v0(:,37),'g')
grid on
title('Voltaje en el lado emisor de la fase C')