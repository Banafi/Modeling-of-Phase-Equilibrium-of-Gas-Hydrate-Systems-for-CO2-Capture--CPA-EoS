function [PHiW]=PHiW_pure(T,P)
global a;global b;global eps;global beta;

R=8.314*(0.01);		%  lit.bar/mol.K
v=1.1*b(1,1); epsilon=0.00000005*v; nnn=0; steplength=10;
while(abs(steplength)>epsilon)
    nnn=nnn+1;
	if(nnn>50)
        nnn;
        break
    end
    poL=1/v;
    eta=b(1,1)*poL/4;
    g=(1-eta/2)/(1-eta)^3;
    delta=g*(exp(eps(1,1)/R/T)-1)*b(1,1)*beta(1,1);
    dlngdv=(2.5-eta)/(1-eta)/(1-eta/2)*(-b(1,1)/4/v^2);
    PPR1=R*T/(v-b(1,1))-a(1,1)/(v*(v+b(1,1))+b(1,1)*(v-b(1,1)));
    XAL=(-1+(1+8*poL*delta)^0.5)/4/poL/delta;       % H2O(4C)
    Pass1=-R*T/2/v*(1-v*dlngdv)*(4-4*XAL);
    Pcal1=PPR1+Pass1;
    Pcal1=P-Pcal1;  
              
	v2=v+0.00000000001*v;
    
    poL=1/v2;
    eta=b(1,1)*poL/4;
    g=(1-eta/2)/(1-eta)^3;
    delta=g*(exp(eps(1,1)/R/T)-1)*b(1,1)*beta(1,1);
    dlngdv=(2.5-eta)/(1-eta)/(1-eta/2)*(-b(1,1)/4/v2^2);
    PPR2=R*T/(v2-b(1,1))-a(1,1)/(v2*(v2+b(1,1))+b(1,1)*(v2-b(1,1)));
    XAL=(-1+(1+8*poL*delta)^0.5)/4/poL/delta;       % H2O(4C)
    Pass2=-R*T/2/v2*(1-v2*dlngdv)*(4-4*XAL);
    Pcal2=PPR2+Pass2;
    Pcal2=P-Pcal2;
         
    dpcaldv=(Pcal2-Pcal1)/(v2-v);
    steplength=Pcal1/dpcaldv;
	vcal=v-0.5*steplength;
    v=vcal;
end
vL=v;

ZL=P*vL/R/T;
poL=1/vL;
eta=b(1,1)*poL/4;
g=(1-eta/2)/(1-eta)^3;
delta=g*(exp(eps(1,1)/R/T)-1)*b(1,1)*beta(1,1);
ARPR=-log(1-b(1,1)/vL)-a(1,1)/(2*2^0.5*b(1,1)*R*T)*log((1+(1+2^0.5)*b(1,1)/vL)/(1+(1-2^0.5)*b(1,1)/vL));
XAL=(-1+(1+8*poL*delta)^0.5)/4/poL/delta;       % H2O(4C)
ARass=4*(log(XAL)-XAL/2)+2;
LnPHiW=ARPR+ARass+ZL-1-log(ZL);
 
PHiW=exp(LnPHiW);
return
