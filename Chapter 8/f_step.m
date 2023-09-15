function [dxi]=step(C,h,u,v,M)
% Function in charge of obtaining an adequate marching step.

Ny=size(u,1);
max=0;
for j=1:Ny
    t1=abs(tan(atan(v(j)/u(j))+(asin(1/M(j)))));
    t2=abs(tan(atan(v(j)/u(j))-(asin(1/M(j)))));
    if(t1>max)
        max=t1;
    end
    if(t2>max)
        max=t2;
    end
end
dxi=C*(h/Ny)/max;
end
