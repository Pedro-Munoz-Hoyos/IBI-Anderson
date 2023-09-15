function [rho] = decoderho(F,gam)
% Function in charge of decoding the density from the flux variables.

A=(0.5*(F(3)^2)/F(1))-F(4);
B=gam*F(1)*F(2)/(gam-1);
C=-0.5*(gam+1)*(F(1)^3)/(gam-1);
rho=(-B+sqrt((B^2)-4*A*C))/(2*A);
end
