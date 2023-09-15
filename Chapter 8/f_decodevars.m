function [u,v,p,T,M] = decodevars(F,rho,gam,R)
% Function in charge of decoding the primitive variables from the flux
% variables.

u=F(1)/rho;
v=F(3)/F(1);
p=F(2)-((F(1)^2)/rho);
T=p/(rho*R);
M=sqrt(((u^2)+(v^2))/(gam*R*T));
end
