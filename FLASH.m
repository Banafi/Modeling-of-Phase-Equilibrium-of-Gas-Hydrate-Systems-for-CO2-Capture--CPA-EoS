function [comp]=FLASH(T,P,z,nc)
Tcp=[647.3 304.25 126.2];    %K
Pcp=[220.60 73.90 33.903];     %bar
omegap=[0.344 0.239 0.039];
for i=1:nc
    knew(i)=exp(log(Pcp(i)/P)+5.373*(1+omegap(i))*(1-Tcp(i)/T));
end
k=[0 0 0];
while abs(knew(1)-k(1))>0.000000001 && abs(knew(2)-k(2))>0.000000001 && abs(knew(3)-k(3))>0.000000001
    k=knew;
    Betamin=0;    Betamax=1;
    j=0;    jj=0;
    for i=1:nc
        if knew(i)>1
            j=j+1;
            A(j)=(knew(i)*z(i)-1)/(knew(i)-1);
            Betamin=max(A);
            if Betamin<0
                Betamin=0;
            else
                Betamin=Betamin;
            end
        else
            jj=jj+1;
            AA(jj)=(1-z(i))/(1-knew(i));
            Betamax=min(AA);
            if Betamax>1
                Betamax=1;
            else
                Betamax=Betamax;
            end
        end
    end
    Beta=0;
    Betanew=(Betamin+Betamax)/2;
    sumy=10;
    sumx=10;
    while (abs(sumy-1)>0.00000001 && abs(sumx-1)>0.00000001)
        Beta=Betanew;
        g=(z(1)*(k(1)-1))/(1+Beta*(k(1)-1))+(z(2)*(k(2)-1))/(1+Beta*(k(2)-1))+(z(3)*(k(3)-1))/(1+Beta*(k(3)-1));
        if g>0
            Betamin=Beta;
        else
            Betamax=Beta;
        end
        Beta=Beta+0.00001;
        gp=(z(1)*(k(1)-1))/(1+Beta*(k(1)-1))+(z(2)*(k(2)-1))/(1+Beta*(k(2)-1))+(z(3)*(k(3)-1))/(1+Beta*(k(3)-1));
        Beta=Beta-0.00001;
        gBeta=(gp-g)/0.00001;
        Betap=Beta-0.5*g/gBeta;
        
        gpp=(z(1)*(k(1)-1))/(1+Betap*(k(1)-1))+(z(2)*(k(2)-1))/(1+Betap*(k(2)-1))+(z(3)*(k(3)-1))/(1+Betap*(k(3)-1));
        
        Betanew=Beta+3*(Betap-Beta)*g/(2*g-gpp);
        if Betanew>Betamin && Betanew<Betamax
            Betanew=Betanew;
        else
            Betanew=(Betamin+Betamax)/2;
        end
        for i=1:nc
        xL(i)=z(i)/(1+Betanew*(k(i)-1));
        xV(i)=k(i)*z(i)/(1+Betanew*(k(i)-1));
        end
        sumx=xL(1)+xL(2)+xL(3);
        sumy=xV(1)+xV(2)+xV(3);
    end
[phiL]=CALPHIL(T,P,xL,nc);
[phiV]=CALPHIV(T,P,xV,nc);
for i=1:nc
    knew(i)=phiL(i)/phiV(i);
end
end
Beta;
g;
[comp]=[xL(1) xL(2) xL(3);xV(1) xV(2) xV(3)];
return