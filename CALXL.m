function [XAL]=CALXL(n,poL,xL,delta)
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%c			Solving  X by iteratiom
%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
for i=1:n
    XAL(i)=0.5;
end
del=10.0;
nn=0;
while(del>0.0001)
		nn=nn+1;
		if(nn>10000)
			nn
			XAL
			break
        end
		XALold=XAL;
%-----------------------------------CO2 as 4C------------------------------;
        XAL(1)=1.0/(1.0+poL*(2*xL(1)*XALold(2)*delta(1,1)+2*xL(2)*XALold(4)*delta(1,2)));
        XAL(2)=1.0/(1.0+poL*(2*xL(1)*XALold(1)*delta(1,1)+2*xL(2)*XALold(3)*delta(2,1)));
        XAL(3)=1.0/(1.0+poL*(2*xL(1)*XALold(2)*delta(1,2)+2*xL(2)*XALold(4)*delta(2,2)));
        XAL(4)=1.0/(1.0+poL*(2*xL(1)*XALold(1)*delta(1,2)+2*xL(2)*XALold(3)*delta(2,2)));
        dif(1)=abs((XAL(1)-XALold(1))/XAL(1));
		dif(2)=abs((XAL(2)-XALold(2))/XAL(2));
		dif(3)=abs((XAL(3)-XALold(3))/XAL(3));
        dif(4)=abs((XAL(4)-XALold(4))/XAL(4));

		del=max(dif);
end
return