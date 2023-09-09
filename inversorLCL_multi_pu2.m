function dx=inversorLCL_multi_pu2(t,x)
%El vector x tiene 13*num+6 variables de estado, donde
%num es la cantidad de inversores en paralelo
%sistema modelado en pu

global Ri Li Rd Cf R L Rf Lf KP KI VCD  MF Lg Rg 
global Ccap kp_pll ki_pll Vp f num
global Vb Sb Ib Zb Anb wb Xb Lb Cb IDREF IQREF



w0=2*pi*f;
vred=Vp*[sin(w0*t);sin(w0*t-2*pi/3);sin(w0*t+2*pi/3)];      
AA=zeros(3,1);
for k=1:num
        ii=x(13*(k-1)+1:13*(k-1)+3);
        vf=x(13*(k-1)+4:13*(k-1)+6);
        ig=x(13*(k-1)+7:13*(k-1)+9);
        gamma=x(13*(k-1)+10);
        wi=x(13*(k-1)+11);
        xind=x(13*(k-1)+12);
        xinq=x(13*(k-1)+13);
        v_cap=x(13*num+1:13*num+3);
        it=x(13*num+4:13*num+6);
        
        Theta1=w0*t+Anb*gamma;

        T1=-2/3*[cos(Theta1) cos(Theta1-2*pi/3) cos(Theta1+2*pi/3)];

        C1=2/3*[sin(Theta1) sin(Theta1-2*pi/3) sin(Theta1+2*pi/3)
                        cos(Theta1) cos(Theta1-2*pi/3) cos(Theta1+2*pi/3)];  

        d_ig1=Rd(k)*wb/(Lf(k)+L(k))*ii-(Rd(k)+Rf(k)+R(k))*wb/(Lf(k)+L(k))*ig+vf*wb/(Lf(k)+L(k))-v_cap*wb/(Lf(k)+L(k));

        vpcc1=Rf(k)*(ig)+Lf(k)/wb*d_ig1+v_cap;

        vtpcc1=T1*vpcc1; %para la sincronizacion1

        %componentes en dq de la corriente y el voltaje
        idq1=C1*ii;
        vred_dq1=C1*vpcc1;
    
        %err_d1=IDREF(k)-idq1(1);
        %err_q1=IQREF(k)-idq1(2);

        %ud1=KP(k)*err_d1/Zb+xind*Xb/Vb;
        %uq1=KP(k)*err_q1/Zb+xinq*Xb/Vb;

        Leq=Li(k)+L(k);

        sd1=2/VCD(k)*(-w0*Leq*idq1(2)/wb+vred_dq1(1)+KP(k)*(IDREF(k)-idq1(1))/Zb+Xb*xind/Vb);
        sq1=2/VCD(k)*(w0*Leq*idq1(1)/wb+vred_dq1(2)+KP(k)*(IQREF(k)-idq1(2))/Zb+Xb*xinq/Vb);

        ma1=sqrt(sd1^2+sq1^2)/Vp;
        theta1=atan2(sq1,sd1);

        %moduladoras para  VSC1
        modA1=ma1*sin(w0*t+theta1);
        modB1=ma1*sin(w0*t+theta1-2*pi/3);
        modC1=ma1*sin(w0*t+theta1+2*pi/3);
        trian1=sawtooth(MF(k)*w0*t+1.5*pi,0.5);
        
        %Aproximacion suave de la SPWM
        sa1=(tanh(99*(modA1-trian1))+1)/2;
        sb1=(tanh(99*(modB1-trian1))+1)/2;
        sc1=(tanh(99*(modC1-trian1))+1)/2;
        ss1=(sa1+sb1+sc1)/3;

        g1=[sa1-ss1;sb1-ss1;sc1-ss1];
        AA=AA+ig;
        dx(13*(k-1)+1:13*k,1)=[-(Ri(k)+Rd(k))*wb/Li(k)*ii-wb/Li(k)*vf+Rd(k)*wb/Li(k)*ig+VCD(k)*wb/Li(k)*g1
                                                wb*1/Cf(k)*ii-wb*1/Cf(k)*ig
                                                d_ig1
                                                -w0/Anb+wi*wb/Anb-kp_pll*vtpcc1*Vb/Anb
                                                 -ki_pll*vtpcc1*Vb/wb
                                                 KI(k)*(IDREF(k)-idq1(1))*Ib/Xb
                                                 KI(k)*(IQREF(k)-idq1(2))*Ib/Xb];                                     
end


dx(13*num+1:13*num+6,1)=[wb*1/(Ccap/Cb)*(AA-it)
                       -(Rg/Zb)*wb/(Lg/Lb)*it+v_cap*wb/(Lg/Lb)-(vred/Vb)*wb/(Lg/Lb)];
                                            
end