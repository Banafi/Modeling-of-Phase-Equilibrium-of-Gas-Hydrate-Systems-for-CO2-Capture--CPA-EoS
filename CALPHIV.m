function [phiV]=CALPHIV(T,P,xV,nc)
global eps; global beta; global Tc; global Pc; global omega; global n; global kij; global b;global a
R=8.314*(0.01);		%  lit.bar/mol.K
for i=1:nc
    m(i)=0.37464+1.54226*omega(i)-0.26992*omega(i)^2;
    C(i)=1+m(i)/2;
    d(i)=(C(i)-1)/C(i);
end
for i=1:nc
    if T/Tc(i)>1
        alfa(i)=(exp(C(i)*(1-(T/Tc(i))^d(i))))^2;
    else
        alfa(i)=(1+m(i)*(1-(T/Tc(i))^0.5))^2;
    end
end
for i=1:nc
    a(i,i)=(0.45724*(R*Tc(i))^2/Pc(i)*alfa(i));
    b(i,i)=(0.0778*R*Tc(i)/Pc(i));
end
for i=1:nc
    for j=1:nc
        a(i,j)=(a(i,i)*a(j,j))^0.5*(1-kij(i,j));      % aij
        b(i,j)=(b(i,i)+b(j,j))/2;       % bij
    end
end
amixV=0;
for i=1:nc
    for j=1:nc
        amixV=amixV+xV(i)*xV(j)*a(i,j);     % a mix. for vapor phase 
    end
end
bmixV=0;
for i=1:nc
    bmixV=bmixV+xV(i)*b(i,i);     % b mix. for vapor phase
end
vV=CALvV(T,P,xV,bmixV,amixV,nc);
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
ZV=P*vV/R/T;
for I=1:nc
    SUMXA=0.0;
	for J=1:nc
		SUMXA=SUMXA+xV(J)*a(I,J);
    end
	
    termprV(I)=-log(P*(vV-bmixV)/R/T) + b(I,I)/bmixV*(P*vV/R/T-1) - (amixV*(-(b(I,I)/bmixV)+(2*(SUMXA))/amixV)*log((1+(1+2^0.5)*bmixV*poV)/(1+(1-2^0.5)*bmixV*poV)))/(2*2^0.5*bmixV*R*T);
%-----------------------------------CO2 as 4C------------------------------;
    TermassocV(I)=b(I,I)/8/g/vV*(2.5-eta)/(1-eta)^4*(2*((1-XAV(3))+(1-XAV(4)))*xV(2)+2*((1-XAV(1))+(1-XAV(2)))*xV(1));
end
%-----------------------------------CO2 as 4C------------------------------;
LNPHI(1)=termprV(1) + TermassocV(1) + 2*(log(XAV(1)) + log(XAV(2))) - log(ZV);
LNPHI(2)=termprV(2) + TermassocV(2) + 2*(log(XAV(3)) + log(XAV(4))) - log(ZV);
LNPHI(3)=termprV(3) + TermassocV(3) - log(ZV);
for i=1:nc
		phiV(i)=exp(LNPHI(i));
end
return