function [T Y]=rk4(f,a,ya,M,u)
%Datos
%       R=rk4(f,a,ya,M)
%       -f es la funci?n almacenada como una cadena de caracteres 'f'
%       -a es el intervalo de soluci?n
%       -ya es la condici?n inicial ya(a)
%       -M es el n?mero de pasos de integraci?n
%       -u es un vector que contiene a las entradas del sistema
% Resultado
%       -R=[T' Y'] siendo T el vector de las abscisas e Y el vector de las
%       ordenadas
M=M-1;
h=(a(2)-a(1))/M;
Y=zeros(length(ya),M+1);
T=a(1):h:a(2);
Y(:,1)=ya;
for j=1:M
    k1=h*feval(f,T(j),Y(:,j),u);
    k2=h*feval(f,T(j)+h/2,Y(:,j)+k1/2,u);
    k3=h*feval(f,T(j)+h/2,Y(:,j)+k2/2,u);
    k4=h*feval(f,T(j)+h,Y(:,j)+k3,u);
    Y(:,j+1)=Y(:,j)+(k1+2*k2+2*k3+k4)/6;
end 
T=T';
Y=Y';
