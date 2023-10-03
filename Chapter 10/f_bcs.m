function [rho,u,v,T,p,Uw,Uo] = bcs(u,v,T,p,rho,cv,R)
% Function in charge of applying adequate boundary conditions in the
% outflow and wall boundaries.

nx=size(u,1);
ny=size(u,1);

% OUTLET
for j=2:ny-1
    u(nx,j)=2*u(nx-1,j)-u(nx-2,j);
    v(nx,j)=2*v(nx-1,j)-v(nx-2,j);
    T(nx,j)=2*T(nx-1,j)-T(nx-2,j);
    p(nx,j)=2*p(nx-1,j)-p(nx-2,j);
    rho(nx,j)=2*rho(nx-1,j)-rho(nx-2,j);
end
Uo=computeU(rho(nx,2:ny-1),cv,u(nx,2:ny-1),v(nx,2:ny-1),T(nx,2:ny-1));

% WALL
for i=2:nx
    p(i,1)=2*p(i,2)-p(i,3);
    %T(i,1)=T(i,2); % UNCOMMENT TO IMPLEMENT AN ADIABATIC WALL B.C.
    rho(i,1)=2*rho(i,2)-rho(i,3);
end
Uw=computeU(rho(2:nx,1),cv,u(2:nx,1),v(2:nx,1),T(2:nx,1));
end

