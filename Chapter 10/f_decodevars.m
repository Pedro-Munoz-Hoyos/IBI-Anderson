function [rho,u,v,T,p] = decodevars(U,cv,R)
% Function in charge of decoding the primitive variables from the flux
% vector U.

nx=size(U,1);
ny=size(U,2);
rho=zeros(nx,ny);
u=zeros(nx,ny);
v=zeros(nx,ny);
T=zeros(nx,ny);
p=zeros(nx,ny);

for i=1:nx
    for j=1:ny
        rho(i,j)=U(i,j,1);
        u(i,j)=U(i,j,2)/U(i,j,1);
        v(i,j)=U(i,j,3)/U(i,j,1);
        T(i,j)=((U(i,j,4)/U(i,j,1))-0.5*((u(i,j)^2)+(v(i,j)^2)))/cv;
        p(i,j)=U(i,j,1)*R*T(i,j);
    end
end

end
