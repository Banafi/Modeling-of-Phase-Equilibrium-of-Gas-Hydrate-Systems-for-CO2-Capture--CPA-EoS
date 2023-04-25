clc, close all, clear all;
format long
%---------for Function VCAL------------------------------------------------;
global eps; global beta; global kij;global Tc;  global Pc;  global omega; global n
R=8.314*(0.01);		%  lit.bar/mol.K
%-----------------------Exp data-------------------------------------------;
zexp1=xlsread('modeltrans.xlsx','D3:D42');     zexp1=zexp1';
zexp2=xlsread('modeltrans.xlsx','E3:E42');     zexp2=zexp2';
Teq=274;
Peq=200;
%Teq=xlsread('CPA-PR3.xlsx','G3:G349');       Teq=Teq';
%Peq=xlsread('CPA-PR3.xlsx','H3:H349');  Peq=Peq*10;      Peq=Peq';
L=length(zexp1);
%--------------------------------------------------------------------------;
fit=[167.95 128.27 2.9663 3.1713];
%fit=[166.909 128.83 2.9665 3.1715];
nc=3;       % number of components
n=4;        % number of associating sites
Tc=[305.40 259.44 124.16];    %K
Pc=[135.62 58.42 30.6];     %bar
omega=[0.1609 0.1290 0.05];
structure=1; langmuirC=2;
%-----------------------CPA Parameters-------------------------------------;
%		   NAME OF COMPONENTS	  
%	  1=water     4C asociation scheme
%	  2=CO2       4C asociation scheme
%	  3=N2        non associating
for nexp=1:L
    nexp
%if zexp2(nexp)==0; structure=2; else structure=1; end
z(1)=zexp1(nexp); z(2)=zexp2(nexp); z(3)=1-z(1)-z(2); T=Teq;
eps=[1811.3*R 0.0 0.0; 0.0 481.1*R 0.0; 0.0 0.0 0.0];
beta=[0.1062 0.0 0.0; 0.0 0.0457 0.0; 0.0 0.0 0.0];
eps(1,2)=(eps(1,1)*eps(2,2))^0.5*(1-0.85180+0.00205*T); eps(2,1)=eps(1,2);
%eps(1,2)=(eps(1,1)+eps(2,2))/2; eps(2,1)=eps(1,2); 
beta(1,2)=(beta(1,1)*beta(2,2))^0.5; beta(2,1)=beta(1,2);
%kij=[0.0 -0.79653+0.00245*T 0.0;-0.79653+0.00245*T 0.0 0.0;0.0 0.0 0.0];
kij=[0.0 -0.79653+0.00245*T 0.37955-350.88/T;-0.79653+0.00245*T 0.0 0.0;0.37955-350.88/T 0.0 0.0];
P=0; 
Pnew=Peq/1;
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
end
%YHYD1=YHYD1'
%YHYD1=xlswrite('CPA-PR3.xlsx',YHYD1,'O347:O349')
PHYD=PHYD'
PHYD=xlswrite('modeltrans.xlsx',PHYD,'O3:O42')


