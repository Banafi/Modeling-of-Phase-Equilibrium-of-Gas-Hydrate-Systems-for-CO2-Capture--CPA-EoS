function HYD=HYDRATE1(T,P,xV,xL,nc,structure,langmuirC)
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%c		       	Parameters for Wan der Waals-Platteeuw model
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
R=8.314;
T0=273.15;      % K
if structure==1
    AAA=4.6e-6;     %DV(B/Ice)+DV(Ice/Liqwater)  (m3/mol)
    BBB=1297;       %Dmio(T0,P0)(B/Liqwater)  (j/mol)
    CCC=1389;      %DH(T0,P0) (B/Ice)  (j/mol)
    DDD=-6011;      %DH(T0,P0) (Ice-Liqwater)  (j/mol)
    aaa=-38.12;
    bbb=-0.1406;
    EEE=aaa-bbb*(T-T0);     %DCp(T)  (j/K/mol)
else
    AAA=5e-6; BBB=937; CCC=1025; DDD=-6011; aaa=-38.12; bbb=-0.141; EEE=aaa-bbb*(T-T0);
end
vm=[1/23 3/23;2/17 1/17];       % vm   vm(i,j) i:structure(1 or 2)  j:cavity type(1=small or 2=large)
%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%c			Calculation of Langmuir adsorption coefficients
%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
if langmuirC==1          % Langmuir adsorption coefficients  C=A/T*exp(B/T)
if structure==1     
    A(1,1)=1.15849e-8; %A2small    K/Pa
    A(1,2)=3.00712e-8; %A3small 
    A(2,1)=0.78920e-7; %A2large 
    A(2,2)=2.18150e-7; %A3large
    B(1,1)=2.86050e3; %B2small   K
    B(1,2)=2.13934e3; %B3small 
    B(2,1)=3.28085e3; %B2large 
    B(2,2)=2.24147e3; %B3large     
else
    A(1,1)=1.26507e-8; %A2small    K/Pa
    A(1,2)=2.83936e-8; %A3small 
    A(2,1)=4.04863e-7; %A2large 
    A(2,2)=7.48338e-7; %A3large
    B(1,1)=2.78974e3; %B2small   K
    B(1,2)=2.17500e3; %B3small 
    B(2,1)=2.82898e3; %B2large 
    B(2,2)=1.86060e3; %B3large
end
for m=1:2       % m:cavity type (1=small or 2=large) 
    for j=1:2       % component (1=CO2 , 2=N2)
        C(m,j)=A(m,j)/T*exp(B(m,j)/T);      %Langmuir adsorption 
    end
end

else              % Langmuir adsorption coefficients  "Kihara"
    
wi=[0.1527533871307258 0.1527533871307258 0.1491729864726037 0.1491729864726037 0.1420961093183820 0.1420961093183820...
0.1316886384491766 0.1316886384491766 0.1181945319615184 0.1181945319615184 0.1019301198172404 0.1019301198172404...
0.0832767415767048 0.0832767415767048 0.0626720483341091 0.0626720483341091 0.0406014298003869 0.0406014298003869...
0.0176140071391521 0.0176140071391521];
xi=[-0.0765265211334973 0.0765265211334973 -0.2277858511416451 0.2277858511416451 -0.3737060887154195 0.3737060887154195...
-0.5108670019508271 0.5108670019508271 -0.6360536807265150 0.6360536807265150 -0.7463319064601508 0.7463319064601508...
-0.8391169718222188 0.8391169718222188 -0.9122344282513259 0.9122344282513259 -0.9639719272779138 0.9639719272779138...
-0.9931285991850949 0.9931285991850949];
kB=1.38064852e-23;
      % (1)=CO2   (2)=N2
akihara(1)=0.6805e-10; akihara(2)=0.3526e-10; epskihara(1)=kB*171.70; epskihara(2)=kB*128.07;
sigmakihara(1)=2.9643e-10; sigmakihara(2)=3.1723e-10;    % ref. herslund thesis 2013
if structure==1  
    Zm(1)=20; Zm(2)=24;   %(1) small cavity, (2) large cavity, Zm: Coordination number
    Rm(1)=3.95e-10; Rm(2)=4.33e-10;  %(1) small cavity, (2) large cavity, Rm: Average cavity radius(m)  
else
    Zm(1)=20; Zm(2)=28;
    Rm(1)=3.91e-10; Rm(2)=4.73e-10;
end
for m=1:2
    for j=1:2
        fint{m,j}=@(r) exp(-(2.*Zm(m).*epskihara(j).*(sigmakihara(j).^12/Rm(m).^11./r.*((1./10.*((1-r./Rm(m)-akihara(j)./...
        Rm(m)).^(-10)-(1+r./Rm(m)-akihara(j)./Rm(m)).^(-10)))+akihara(j)./Rm(m).*(1./11.*((1-r./Rm(m)-akihara(j)./...
        Rm(m)).^(-11)-(1+r./Rm(m)-akihara(j)./Rm(m)).^(-11))))-sigmakihara(j).^6./Rm(m).^5./r.*((1./4.*((1-r./Rm(m)...
        -akihara(j)./Rm(m)).^(-4)-(1+r./Rm(m)-akihara(j)./Rm(m)).^(-4)))+akihara(j)./Rm(m).*(1./5.*((1-r./Rm(m)-...
        akihara(j)./Rm(m)).^(-5)-(1+r./Rm(m)-akihara(j)./Rm(m)).^(-5))))))./kB./T).*r.^2;
        bint(m,j)=Rm(m)-akihara(j);
        aint=0;
    end
end
for m=1:2
    for j=1:2
        int=0;
        for i=1:20
            int=int+wi(i)*fint{m,j}((bint(m,j)-aint)/2*xi(i)+(aint+bint(m,j))/2);
        end
        int1=(bint(m,j)-aint)/2*int;
        C(m,j)=4*pi/kB/T*int1;
    end
end
end
%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%c			Calculation of Pressure with Wan der Waals-Platteeuw model
%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Pnew1=P*1e5;
P=0;
epsilon1=0.000000000001*Pnew1;
steplength1=10;
while abs(P-Pnew1)>epsilon1
    P=Pnew1;
    P=P*1e-5;
    [phiV]=CALPHIV(T,P,xV,nc);
    phiwpure=PHiW_pure(T,P);
    [phiL]=CALPHIL(T,P,xL,nc);
    P=P*1e5;
    for j=1:2       % fugacity of guest components(1=CO2 , 2=N2)
        fug(j)=phiV(j+1)*xV(j+1)*P;
    end
    sum=0;
    if structure==1
        for m=1:2
            sum1=0;
            for j=1:2
                sum1=sum1+C(m,j)*fug(j);
            end
            sum=sum+vm(1,m)*log(1+sum1);
        end
    else
        for m=1:2
            sum1=0;
            for j=1:2
                sum1=sum1+C(m,j)*fug(j);
            end
            sum=sum+vm(2,m)*log(1+sum1);
        end
    end
    %mio1=sum;
    mio1=sum+log(xL(1)*phiL(1)/phiwpure);
    mio2=(BBB/R/T0)-((CCC+DDD-aaa*T0-bbb*T0^2)/R*(1/T0-1/T)+(aaa+2*bbb*T0)/R*log(T/T0)-bbb/R*(T-T0))+(AAA/R/T*P);
    dif1=mio2-mio1;
    
    P2=P+0.00000001*P;
    P2=P2*1e-5;
    [phiV]=CALPHIV(T,P2,xV,nc);
    phiwpure=PHiW_pure(T,P2);
    [phiL]=CALPHIL(T,P2,xL,nc);
    P2=P2*1e5;
    for j=1:2       % fugacity of guest components(1=CO2 , 2=N2)
        fug(j)=phiV(j+1)*xV(j+1)*P2;
    end
    sum=0;
    if structure==1
        for m=1:2
            sum1=0;
            for j=1:2
                sum1=sum1+C(m,j)*fug(j);
            end
            sum=sum+vm(1,m)*log(1+sum1);
        end
    else
        for m=1:2
            sum1=0;
            for j=1:2
                sum1=sum1+C(m,j)*fug(j);
            end
            sum=sum+vm(2,m)*log(1+sum1);
        end
    end
    %mio1p=sum;
    mio1p=sum+log(xL(1)*phiL(1)/phiwpure);
    mio2p=(BBB/R/T0)-((CCC+DDD-aaa*T0-bbb*T0^2)/R*(1/T0-1/T)+(aaa+2*bbb*T0)/R*log(T/T0)-bbb/R*(T-T0))+(AAA/R/T*P2);
    dif2=mio2p-mio1p;
    
    steplength1=dif1/((dif2-dif1)/(P2-P));
    Pnew=P-0.5*steplength1;
    
    Pnew=Pnew*1e-5;
    [phiV]=CALPHIV(T,Pnew,xV,nc);
    phiwpure=PHiW_pure(T,Pnew);
    [phiL]=CALPHIL(T,Pnew,xL,nc);
    Pnew=Pnew*1e5;
    for j=1:2       % fugacity of guest components(1=CO2 , 2=N2)
        fug(j)=phiV(j+1)*xV(j+1)*Pnew;
    end
    sum=0;
    if structure==1
        for m=1:2
            sum1=0;
            for j=1:2
                sum1=sum1+C(m,j)*fug(j);
            end
            sum=sum+vm(1,m)*log(1+sum1);
        end
    else
        for m=1:2
            sum1=0;
            for j=1:2
                sum1=sum1+C(m,j)*fug(j);
            end
            sum=sum+vm(2,m)*log(1+sum1);
        end
    end
    %mio1pp=sum;
    mio1pp=sum+log(xL(1)*phiL(1)/phiwpure);
    mio2pp=(BBB/R/T0)-((CCC+DDD-aaa*T0-bbb*T0^2)/R*(1/T0-1/T)+(aaa+2*bbb*T0)/R*log(T/T0)-bbb/R*(T-T0))+(AAA/R/T*Pnew);
    difnew=mio2pp-mio1pp;
    
    Pnew1=P+3*(Pnew-P)*dif1/(2*dif1-difnew);
    
end
 difnew;
%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%c			Calculation of water free hydrate composition
%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
for m=1:2       % m:cavity type (1=small or 2=large)
    for j=1:2   % component (1=CO2 , 2=N2)
        sum2=0;
        for L=1:2
            sum2=sum2+C(m,L)*fug(L);
        end
        teta(m,j)=C(m,j)*fug(j)/(1+sum2);
    end
end
if structure==1
    sum3=0;
    for m=1:2
        sum4=0;
        for j=1:2
            sum4=sum4+teta(m,j);
        end
        sum3=sum3+vm(1,m)*sum4;
    end
else
    sum3=0;
    for m=1:2
        sum4=0;
        for j=1:2
            sum4=sum4+teta(m,j);
        end
        sum3=sum3+vm(2,m)*sum4;
    end
end
if structure==1
    for j=1:2
        sum5=0;
        for m=1:2
            sum5=sum5+vm(1,m)*teta(m,j)/sum3;
        end
        yHYD(j)=sum5;
    end
else
    for j=1:2
        sum5=0;
        for m=1:2
            sum5=sum5+vm(2,m)*teta(m,j)/sum3;
        end
        yHYD(j)=sum5;
    end
end             
PHYD=Pnew*1e-5;
HYD=[PHYD yHYD(1) yHYD(2)];
return