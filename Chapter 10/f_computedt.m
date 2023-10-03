function [dt] = computedt(rho,mu,u,v,T,gam,R,dx,dy,K)
% Function in charge of obtaining an adequate time-step based on the
% following version of the CFL criterion (pag.457).

nx=size(u,1);
ny=size(u,2);

maxvp=0;
for i=2:nx-1
    for j=2:ny-1
        vp=4*gam*(mu(i,j)^2)/(3*0.71*rho(i,j));
        if vp>maxvp
            maxvp=vp;
        end
    end
end

mindt=1000;
for i=2:nx-1
    for j=2:ny-1
        A=(abs(u(i,j))/dx)+(abs(v(i,j))/dy);
        B=sqrt(gam*R*T(i,j))*sqrt(((1/dx)^2)+((1/dy)^2));
        C=2*maxvp*(((1/dx)^2)+((1/dy)^2));
        dt=1/(A+B+C);
        if dt<mindt
            mindt=dt;
        end
    end
end
dt=mindt*K;

end
