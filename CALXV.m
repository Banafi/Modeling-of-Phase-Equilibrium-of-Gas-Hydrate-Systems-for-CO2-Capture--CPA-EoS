function [XAV]=CALXV(n,poV,xV,delta)
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%c			Solving  X by iteratiom
%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
for i=1:n
    XAV(i)=0.5;
end
del=10.0;
nn=0;
while(del>0.0001)
		nn=nn+1;
		if(nn>10000)
			nn
			XAV
			break
        end
		XAVold=XAV;
%-----------------------------------CO2 as 4C------------------------------;
        XAV(1)=1.0/(1.0+poV*(2*xV(1)*XAVold(2)*delta(1,1)+2*xV(2)*XAVold(4)*delta(1,2)));
        XAV(2)=1.0/(1.0+poV*(2*xV(1)*XAVold(1)*delta(1,1)+2*xV(2)*XAVold(3)*delta(2,1)));
        XAV(3)=1.0/(1.0+poV*(2*xV(1)*XAVold(2)*delta(1,2)+2*xV(2)*XAVold(4)*delta(2,2)));
        XAV(4)=1.0/(1.0+poV*(2*xV(1)*XAVold(1)*delta(1,2)+2*xV(2)*XAVold(3)*delta(2,2)));
        dif(1)=abs((XAV(1)-XAVold(1))/XAV(1));
		dif(2)=abs((XAV(2)-XAVold(2))/XAV(2));
		dif(3)=abs((XAV(3)-XAVold(3))/XAV(3));
        dif(4)=abs((XAV(4)-XAVold(4))/XAV(4));


		del=max(dif);
end
return