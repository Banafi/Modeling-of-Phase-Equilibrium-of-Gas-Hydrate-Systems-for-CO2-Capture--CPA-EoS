function [vL]=CALvL(T,P,xL,bmixL,amixL,nc)
global eps; global beta; global b; global n;
%----------------------Parameters for CO2 as 4C----------------------------;
R=8.314*(0.01);		%  lit.bar/mol.K
vLcalnew=1.1*bmixL; vL=0;
epsilon=0.00000005*vLcalnew;
nnn=0; steplength=10;
while(abs(vL-vLcalnew)>epsilon)
    vL=vLcalnew;
    nnn=nnn+1;
	if(nnn>100)
        nnn
        break
    end
    poL=1./vL;
	eta=bmixL*poL/4.;
	g=(1-eta/2)/(1-eta)^3;
    %g=1./(1.-1.9*eta);
    for i=1:nc					
        for j=1:nc
            delta(i,j)=(exp(eps(i,j)/(R*T))-1)*g*beta(i,j)*b(i,j);
        end
    end
    [XAL]=CALXL(n,poL,xL,delta);  % Call solver
    dlngdv=(2.5-eta)/(1-eta)/(1-eta/2)*(-bmixL/4/vL^2);
    %dlngdv=(-0.475*bmixL)/((1.-(0.475*bmixL)/vL)*vL^2);
%----------------------------------CO2 as 4--------------------------------;
    pcal1=(R*T)/(-bmixL+vL)-amixL/(vL*(bmixL+vL)+bmixL*(vL-bmixL))-(0.5*R*T*(1.- dlngdv*vL)*(2*xL(1)*(1-XAL(1)+1-XAL(2))+2*xL(2)*(1-XAL(3)+1-XAL(4))))/vL;
    
    pcal1=P-pcal1;
    
	vL2=vL+0.00000001*vL;
    
    poL=1./vL2;
	eta=bmixL*poL/4.;
	g=(1-eta/2)/(1-eta)^3;
    %g=1./(1.-1.9*eta);
    for i=1:nc					
        for j=1:nc
            delta(i,j)=(exp(eps(i,j)/(R*T))-1)*g*beta(i,j)*b(i,j);
        end
    end
    [XAL]=CALXL(n,poL,xL,delta);  % Call solver
    dlngdv=(2.5-eta)/(1-eta)/(1-eta/2)*(-bmixL/4/vL2^2);
    %dlngdv=(-0.475*bmixL)/((1.-(0.475*bmixL)/vL2)*vL2^2);
%-----------------------------CO2 as 4C------------------------------------;
    pcal2=(R*T)/(-bmixL+vL2)-amixL/(vL2*(bmixL+vL2)+bmixL*(vL2-bmixL))-(0.5*R*T*(1.- dlngdv*vL2)*(2*xL(1)*(1-XAL(1)+1-XAL(2))+2*xL(2)*(1-XAL(3)+1-XAL(4))))/vL2;
    pcal2=P-pcal2;
    
    dpcaldvL=(pcal2-pcal1)/(vL2-vL);
    steplength=pcal1/dpcaldvL;
	vLcal=vL-0.5*steplength;
    
    poL=1./vLcal;
	eta=bmixL*poL/4.;
	g=(1-eta/2)/(1-eta)^3;
    %g=1./(1.-1.9*eta);
    for i=1:nc					
        for j=1:nc
            delta(i,j)=(exp(eps(i,j)/(R*T))-1)*g*beta(i,j)*b(i,j);
        end
    end
    [XAL]=CALXL(n,poL,xL,delta);  % Call solver
    dlngdv=(2.5-eta)/(1-eta)/(1-eta/2)*(-bmixL/4/vLcal^2);
    %dlngdv=(-0.475*bmixL)/((1.-(0.475*bmixL)/vL2)*vL2^2);
%-----------------------------CO2 as 4C------------------------------------;
    pcalnew=(R*T)/(-bmixL+vLcal)-amixL/(vLcal*(bmixL+vLcal)+bmixL*(vLcal-bmixL))-(0.5*R*T*(1.- dlngdv*vLcal)*(2*xL(1)*(1-XAL(1)+1-XAL(2))+2*xL(2)*(1-XAL(3)+1-XAL(4))))/vLcal;
    pcalnew=P-pcalnew;
    
    vLcalnew=vL+3*(vLcal-vL)*pcal1/(2*pcal1-pcalnew);  
end
pcalnew;
return