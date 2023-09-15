function [S] = artivis(Fl,Fc,Fr,rho,Cy)
% Function in charge of obtaining the artificial viscosity.

S=zeros(4,1);
spl=Fl(2)-((Fl(1)^2)/rho);
sp=Fc(2)-((Fc(1)^2)/rho);
spr=Fr(2)-((Fr(1)^2)/rho);
for i=1:4
    S(i)=Cy*abs(spl-2*sp+spr)*(Fl(i)-2*Fc(i)+Fr(i))/(spl+2*sp+spr);
end
end
