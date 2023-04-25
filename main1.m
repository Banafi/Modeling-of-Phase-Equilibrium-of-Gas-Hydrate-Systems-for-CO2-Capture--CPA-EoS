clc, close all, clear all;
format long
%---------for Function VCAL------------------------------------------------;
global eps; global beta; global kij; global Tc;  global Pc;  global omega; global n
R=8.314*(0.01);		%  lit.bar/mol.K
%-----------------------Exp data-------------------------------------------;
%zexp1=xlsread('CPA-PR3.xlsx','D4:D57');     zexp1=zexp1';
%zexp2=xlsread('CPA-PR3.xlsx','E4:E57');     zexp2=zexp2';
%Teq=xlsread('CPA-PR3.xlsx','G4:G57');       Teq=Teq';
%Peq=xlsread('CPA-PR3.xlsx','H4:H57');  Peq=Peq*10;      Peq=Peq';
%L=length(Peq);
T=273;
y1=0.06;
z(1)=0.8;
z(2)=(1-z(1))*y1;
z(3)=1-z(1)-z(2);
%--------------------------------------------------------------------------;
%nexp=53;
fit=[167.95 128.27 2.9663 3.1713];
nc=3;       % number of components
n=4;        % number of associating sites
structure=2;
langmuirC=2;
%z(1)=zexp1(nexp);
%z(2)=zexp2(nexp);
%z(3)=1-z(1)-z(2);
%T=Teq(nexp);
%-----------------------CPA Parameters-------------------------------------;
%		   NAME OF COMPONENTS	  
%	  1=water     4C asociation scheme
%	  2=CO2       4C asociation scheme
%	  3=N2        non associating
eps=[1811.3*R 0.0 0.0; 0.0 481.1*R 0.0; 0.0 0.0 0.0];
beta=[0.1062 0.0 0.0; 0.0 0.0457 0.0; 0.0 0.0 0.0];
eps(1,2)=(eps(1,1)*eps(2,2))^0.5*(1-0.85180+0.00205*T);
%eps(1,2)=(eps(1,1)+eps(2,2))/2;
eps(2,1)=eps(1,2);
beta(1,2)=(beta(1,1)*beta(2,2))^0.5;
beta(2,1)=beta(1,2);
%-----------------------Tc, Pc, omega Parameters---------------------------;
Tc=[305.40 259.44 124.16];    %K
Pc=[135.62 58.42 30.6];     %bar
omega=[0.1609 0.1290 0.05];
%kij=[0.0 -0.79653+0.00245*T 0.0;-0.79653+0.00245*T 0.0 0.0;0.0 0.0 0.0];
kij=[0.0 -0.79653+0.00245*T 0.37955-350.88/T;-0.79653+0.00245*T 0.0 0.0;0.37955-350.88/T 0.0 0.0];
P=0;
Pnew=300;
%Pnew=Peq(nexp)/2;       %bar
while (abs(Pnew-P)>0.000001)
P=Pnew;
[comp]=FLASH(T,P,z,nc);
xL=comp(1,:);
xV=comp(2,:);
HYD=HYDRATE(fit,T,P,xV,xL,nc,structure,langmuirC);
Pnew=HYD(1)
end
Xeq=xL
Yeq=xV
YHYD=[HYD(2) HYD(3)]
PHYD=Pnew


