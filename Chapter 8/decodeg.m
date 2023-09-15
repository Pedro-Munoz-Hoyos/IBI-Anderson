function [G] = decodeg(F,rho,gam)
% Function in charge of decoding the flux vector G from the flux vector F.

G=zeros(4,1);
G(1)=rho*F(3)/F(1);
G(2)=F(3);
G(3)=(rho*((F(3)/F(1))^2))+F(2)-((F(1)^2)/rho);
G(4)=((gam*(F(2)-((F(1)^2)/rho))/(gam-1))+0.5*rho*(((F(1)/rho)^2)+((F(3)/F(1))^2)))*F(3)/F(1);
end
