function [phiL]=CALPHIL(T,P,xL,nc)
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
amixL=0;
for i=1:nc
    for j=1:nc
        amixL=amixL+xL(i)*xL(j)*a(i,j);     % a mix. for liquid phase 
    end
end
bmixL=0;
for i=1:nc
    bmixL=bmixL+xL(i)*b(i,i);     % b mix. for liquid phase
end
vL=CALvL(T,P,xL,bmixL,amixL,nc);
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
ZL=P*vL/R/T;
for I=1:nc
    SUMXA=0.0;
	for J=1:nc
		SUMXA=SUMXA+xL(J)*a(I,J);
    end
	termprL(I)=-log(P*(vL-bmixL)/R/T) + b(I,I)/bmixL*(P*vL/R/T-1) - (amixL*(-(b(I,I)/bmixL)+(2*(SUMXA))/amixL)*log((1+(1+2^0.5)*bmixL*poL)/(1+(1-2^0.5)*bmixL*poL)))/(2*2^0.5*bmixL*R*T);
%-----------------------------------CO2 as 4C------------------------------;
    TermassocL(I)=b(I,I)/8/g/vL*(2.5-eta)/(1-eta)^4*(2*((1-XAL(3))+(1-XAL(4)))*xL(2)+2*((1-XAL(1))+(1-XAL(2)))*xL(1));
end
%-----------------------------------CO2 as 4C------------------------------;
LNPHI(1)=termprL(1) + TermassocL(1) + 2*(log(XAL(1)) + log(XAL(2))) - log(ZL);
LNPHI(2)=termprL(2) + TermassocL(2) + 2*(log(XAL(3)) + log(XAL(4))) - log(ZL);
LNPHI(3)=termprL(3) + TermassocL(3) - log(ZL);
for i=1:nc
		phiL(i)=exp(LNPHI(i));
end
return

