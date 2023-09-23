function [Z,Y]=TL_param_Nhaz_JULIO(s,Pos,xx)
  %Calculo de las matrices de impedancia y admitancia para una linea que tiene
  %conductores en paralelo por fase ()semejantes a un haz
  % xx es un vector que tiene los siguientes datos (en el orden mencionado):
  % 1.-Radio del conductor (en metros) 
  % 2.-Resistencia de cd en Ohm/m
  % 3.-Numero de conductores totales
  % 4.-Resistividad del terreno en Ohms*m
  % 5.-Hay arreglo de haz simetrico por fase? (S=1/N=0)
  % 6.-Cantidad de conductores por haz (Si no hay poner 1)
  % 7.-Diametro del circulo equivalente sobre el cual se distribuye el haz en m(Si no hay poner 0)
  
  mu0=4*pi*1e-7;
  eps0=8.854e-12;
  r_cond=xx(1);
  Rcd=xx(2);
  Ncond=xx(3);
  rho_earth=xx(4);
  nn=xx(6);
  R=xx(7)/2;
  rho_cond=Rcd*pi*r_cond^2;
  
%Constantes relativas de algo, no tocar (es como la mu relativa del material)
  murc = 1;
  murt = 1;
  
  x=real(Pos)';
  y=imag(Pos)';
  
pcom = sqrt(1/(s*mu0*murt/rho_earth)); %complex penetration depth of ground plane  
r_cond=(nn*R^(nn-1)*r_cond)^(1/nn); %Redefinicion del radio equivalente
%Calculation of GMR and GMD for general bundle case (no symmetry assumed)
DI=zeros(Ncond);%Numerador
DR=zeros(Ncond);%Denominador
DIp=zeros(Ncond);%Numerador p
DRp=zeros(Ncond);%Denominador p
for k1=1:Ncond
  for k2=1:Ncond
    if k1==k2
      DI(k1,k2)=2*y(k1);
      DR(k1,k2)=r_cond;
      
      DIp(k1,k2)=2*(y(k1)+pcom)-r_cond;
      DRp(k1,k2)=2*y(k1)-r_cond;
    else
      DI(k1,k2)=sqrt((y(k1)+y(k2))^2+(x(k1)-x(k2))^2);
      DR(k1,k2)=sqrt((y(k1)-y(k2))^2+(x(k1)-x(k2))^2);
      
      DIp(k1,k2)=sqrt((y(k1)+y(k2)+2*pcom)^2+(x(k1)-x(k2))^2);
      DRp(k1,k2)=sqrt((y(k1)+y(k2))^2+(x(k1)-x(k2))^2);
    end
  end
end

r_cond=xx(1);
POT = log( DI./DR);
Zaf = (0.5/pi)*sqrt(s*mu0*murc*rho_cond)/r_cond;
Zcond = sqrt(Rcd^2+Zaf^2)/nn;

%Series impedance and shunt admittance matrices per unit length
Z = (s*mu0/(2*pi))*POT+(Zcond)*eye(Ncond)+(s*mu0*murt/(2*pi))*log(DIp./DRp);
Y = s*2*pi*eps0*inv(POT);  
  
end
