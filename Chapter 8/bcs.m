function [rhoact,uact,vact,pact,Tact,Mact,F,G] = bcs(u,v,p,T,M,gam,R,wang)
% Function in charge of applying the appropiate boundary condition
% according to Abbett's method.

phi=wang-atan(abs(v/u));
fcal=sqrt((gam+1)/(gam-1))*atan(sqrt((gam-1)*((M^2)-1)/(gam+1)))-atan(sqrt((M^2)-1));
fact=fcal+phi;

Mact=1;
err=1;
while err>10^-7
    err=abs(sqrt((gam+1)/(gam-1))*atan(sqrt((gam-1)*((Mact^2)-1)/(gam+1)))-atan(sqrt((Mact^2)-1))-fact);
    Mact=Mact+10^-7;
end

pact=p*(((1+0.5*(gam-1)*(M^2))/(1+0.5*(gam-1)*(Mact^2)))^(gam/(gam-1)));
Tact=T*(1+0.5*(gam-1)*(M^2))/(1+0.5*(gam-1)*(Mact^2));
rhoact=pact/(R*Tact);
uact=Mact*sqrt(gam*R*Tact)*cos(wang);
vact=-uact*tan(wang);

F(1)=rhoact*uact;
F(2)=F(1)*uact+pact;
F(3)=rhoact*uact*vact;
F(4)=(gam*pact*uact/(gam-1))+F(1)*0.5*((uact^2)+(vact^2));
[G]=decodeg(F,rhoact,gam);
end
