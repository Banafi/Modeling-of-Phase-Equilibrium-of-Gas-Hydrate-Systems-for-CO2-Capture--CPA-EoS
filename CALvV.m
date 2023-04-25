function [vV]=CALvV(T,P,xV,bmixV,amixV,nc)
global eps; global beta; global b; global n;
%----------------------Parameters for CO2 as 4C----------------------------;
R=8.314*(0.01);		%  lit.bar/mol.K
%vVcalnew=1.1*bmixV;
vVcalnew=R*T/P; vV=0;
epsilon=0.0000005*vVcalnew;
nnnvv=0;
steplength=10;
while(abs(vV-vVcalnew)>epsilon)
    vV=vVcalnew;
    nnnvv=nnnvv+1;
	if(nnnvv>100)
        nnnvv
        break
    end
    poV=1./vV;
	eta=bmixV*poV/4.;
	g=(1-eta/2)/(1-eta)^3;
    %g=1./(1.-1.9*eta);
    for i=1:nc					
        for j=1:nc
            delta(i,j)=(exp(eps(i,j)/(R*T))-1)*g*beta(i,j)*b(i,j);
        end
    end
    [XAV]=CALXV(n,poV,xV,delta);  % Call solver
    dlngdv=(2.5-eta)/(1-eta)/(1-eta/2)*(-bmixV/4/vV^2);
    %dlngdv=(-0.475*bmixV)/((1.-(0.475*bmixV)/vV)*vV^2);
%-----------------------------CO2 as 4C------------------------------------;
    pcal1=(R*T)/(-bmixV+vV)-amixV/(vV*(bmixV+vV)+bmixV*(vV-bmixV))-(0.5*R*T*(1.- dlngdv*vV)*(2*xV(1)*(1-XAV(1)+1-XAV(2))+2*xV(2)*(1-XAV(3)+1-XAV(4))))/vV;
    pcal1=P-pcal1;
    
	vV2=vV+0.00000001*vV;
    
    poV=1./vV2;
	eta=bmixV*poV/4.;
	g=(1-eta/2)/(1-eta)^3;
    %g=1./(1.-1.9*eta);
    for i=1:nc					
        for j=1:nc
            delta(i,j)=(exp(eps(i,j)/(R*T))-1)*g*beta(i,j)*b(i,j);
        end
    end
    [XAV]=CALXV(n,poV,xV,delta);  % Call solver
    dlngdv=(2.5-eta)/(1-eta)/(1-eta/2)*(-bmixV/4/vV2^2);
    %dlngdv=(-0.475*bmixV)/((1.-(0.475*bmixV)/vV2)*vV2^2);
%-----------------------------CO2 as 4C------------------------------------;
    pcal2=(R*T)/(-bmixV+vV2)-amixV/(vV2*(bmixV+vV2)+bmixV*(vV2-bmixV))-(0.5*R*T*(1.- dlngdv*vV2)*(2*xV(1)*(1-XAV(1)+1-XAV(2))+2*xV(2)*(1-XAV(3)+1-XAV(4))))/vV2;
    pcal2=P-pcal2;
    
    dpcaldvV=(pcal2-pcal1)/(vV2-vV);
    steplength=pcal1/dpcaldvV;
	vVcal=vV-0.5*steplength;
    
    poV=1./vVcal;
	eta=bmixV*poV/4.;
	g=(1-eta/2)/(1-eta)^3;
    %g=1./(1.-1.9*eta);
    for i=1:nc					
        for j=1:nc
            delta(i,j)=(exp(eps(i,j)/(R*T))-1)*g*beta(i,j)*b(i,j);
        end
    end
    [XAV]=CALXV(n,poV,xV,delta);  % Call solver
    dlngdv=(2.5-eta)/(1-eta)/(1-eta/2)*(-bmixV/4/vVcal^2);
    %dlngdv=(-0.475*bmixV)/((1.-(0.475*bmixV)/vV2)*vV2^2);
%-----------------------------CO2 as 4C------------------------------------;
    pcalnew=(R*T)/(-bmixV+vVcal)-amixV/(vVcal*(bmixV+vVcal)+bmixV*(vVcal-bmixV))-(0.5*R*T*(1.- dlngdv*vVcal)*(2*xV(1)*(1-XAV(1)+1-XAV(2))+2*xV(2)*(1-XAV(3)+1-XAV(4))))/vVcal;
    pcalnew=P-pcalnew;
    
    vVcalnew=vV+3*(vVcal-vV)*pcal1/(2*pcal1-pcalnew);
end
%pcalnewVVV=pcalnew;
return