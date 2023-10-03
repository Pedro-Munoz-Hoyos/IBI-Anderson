function [mu] = computemu(mu_fs,T_fs,T)
% Function in charge of obtainin the dynamic viscosity according to
% Sutherland's law.

nx=size(T,1);
ny=size(T,2);
mu=zeros(nx,ny);

for i=1:nx
    for j=1:ny
        mu(i,j)=mu_fs*((T(i,j)/T_fs)^1.5)*(T_fs+110)/(T(i,j)+110);
    end
end


end
