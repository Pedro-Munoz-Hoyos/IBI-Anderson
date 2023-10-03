function [U] = computeU(rho,cv,u,v,T)
% Function in charge of obtaining the flux vector U from the primitive
% variables.

nx=size(u,1);
ny=size(u,2);
U=zeros(nx,ny,4);

for i=1:nx
    for j=1:ny
        U(i,j,1)=rho(i,j);
        U(i,j,2)=rho(i,j)*u(i,j);
        U(i,j,3)=rho(i,j)*v(i,j);
        U(i,j,4)=rho(i,j)*((cv*T(i,j))+0.5*((u(i,j)^2)+(v(i,j)^2)));
    end
end

end
