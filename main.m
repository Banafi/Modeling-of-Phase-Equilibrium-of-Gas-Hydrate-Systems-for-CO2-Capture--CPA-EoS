clc, close all, clear all;
format long
%---------for Function VCAL------------------------------------------------;
global eps; global beta; global kij;global Tc;  global Pc;  global omega; global n
R=8.314*(0.01);		%  lit.bar/mol.K
%-----------------------Exp data-------------------------------------------;
zexp1=xlsread('CPA-PR3.xlsx','D252:D261');     zexp1=zexp1';
zexp2=xlsread('CPA-PR3.xlsx','E252:E261');     zexp2=zexp2';
Teq=xlsread('CPA-PR3.xlsx','G252:G261');       Teq=Teq';
Peq=xlsread('CPA-PR3.xlsx','H252:H261');  Peq=Peq*10;      Peq=Peq';
L=length(Peq);
%--------------------------------------------------------------------------;
%fit=[167.95 128.27 2.9663 3.1713];
%fit=[172.515172 128.5875 3.286787 3.0375];   % fit ba pure
%fit=[167.95 128.5875 2.96637 3.0375];     % CO2 asli , N2 fit ba pure
%fit=[172.515172 128.27 3.28678 3.1713];    % N2 asli , CO2 fit ba pure
fit=[172.315172 128.5875 3.290787 3.0375];
nc=3;       % number of components
n=4;        % number of associating sites
Tc=[305.40 259.44 124.16];    %K
Pc=[135.62 58.42 30.6];     %bar
omega=[0.1609 0.1290 0.05];
langmuirC=2;
%-----------------------CPA Parameters-------------------------------------;
%		   NAME OF COMPONENTS	  
%	  1=water     4C asociation scheme
%	  2=CO2       4C asociation scheme
%	  3=N2        non associating
for nexp=1:L
if zexp2(nexp)==0; structure=2; else structure=1; end
nexp
z(1)=zexp1(nexp); z(2)=zexp2(nexp); z(3)=1-z(1)-z(2); T=Teq(nexp);
eps=[1811.3*R 0.0 0.0; 0.0 481.1*R 0.0; 0.0 0.0 0.0];
beta=[0.1062 0.0 0.0; 0.0 0.0457 0.0; 0.0 0.0 0.0];
eps(1,2)=(eps(1,1)*eps(2,2))^0.5*(1-0.85180+0.00205*T); eps(2,1)=eps(1,2);
%eps(1,2)=(eps(1,1)+eps(2,2))/2; eps(2,1)=eps(1,2); 
beta(1,2)=(beta(1,1)*beta(2,2))^0.5; beta(2,1)=beta(1,2);
%kij=[0.0 -0.79653+0.00245*T 0.0;-0.79653+0.00245*T 0.0 0.0;0.0 0.0 0.0];
kij=[0.0 -0.79653+0.00245*T 0.37955-350.88/T;-0.79653+0.00245*T 0.0 0.0;0.37955-350.88/T 0.0 0.0];
P=0; 
%if nexp==28; Pnew=Peq(nexp)/2;   else  Pnew=Peq(nexp)/1; end     %Kim,2011
%if nexp>33; Pnew=Peq(nexp)/2; else Pnew=Peq(nexp); end         %Kang 2001-1
%if nexp==6; Pnew=Peq(nexp)/2; else Pnew=Peq(nexp); end         %Seo 2002-2
%Pnew=Peq(nexp)/2;                                              %Fan&Guo 1991-1
%if nexp<24; Pnew=Peq(nexp)/2;   elseif nexp<27  Pnew=Peq(nexp)/1; else Pnew=Peq(nexp)/0.5; end
%if nexp<18; Pnew=Peq(nexp)/1; else Pnew=Peq(nexp)/2; end   % kim 2011
Pnew=Peq(nexp)/2;
while (abs(Pnew-P)>0.000001)
P=Pnew;
[comp]=FLASH(T,P,z,nc);
xL=comp(1,:); xV=comp(2,:);
HYD=HYDRATE(fit,T,P,xV,xL,nc,structure,langmuirC);
Pnew=HYD(1);
end
Xeq=xL;
Yeq=xV;
YHYD=[HYD(2) HYD(3)];
PHYD(nexp)=Pnew;
YHYD1(nexp)=YHYD(1);
Yeq1(nexp)=Yeq(2);
SF(nexp)=(YHYD(1)/YHYD(2))/(Yeq(2)/Yeq(3));
end
%SF=SF'
%SF=xlswrite('CPA-PR3.xlsx',SF,'T351:T361')
%Yeq1=Yeq1'
%Yeq1=xlswrite('CPA-PR3.xlsx',Yeq1,'Q362:Q399')
%YHYD1=YHYD1'
%YHYD1=xlswrite('CPA-PR3.xlsx',YHYD1,'S351:S387')
PHYD=PHYD'
PHYD=xlswrite('CPA-PR3.xlsx',PHYD,'Q252:Q261')



