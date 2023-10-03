function [F] = computeF(rho,mu,cv,cp,dx,dy,u,v,T,p,pred)
% Function in charge of obtaining the flux vector F from the primitive
% variables, taking into account whether it is being calculated in the 
% predictor or corrector step in order to mantain second-order accuracy.

nx=size(u,1);
ny=size(u,2);
F=zeros(nx,ny,4);
[tau,q]=tau_q(mu,cp,u,v,T,dx,dy,pred,0);

for i=1:nx
    for j=1:ny
        F(i,j,1)=rho(i,j)*v(i,j);
        F(i,j,2)=F(i,j,1)*u(i,j)-tau(i,j,1);
        F(i,j,3)=F(i,j,1)*v(i,j)+p(i,j)-tau(i,j,2);
        F(i,j,4)=(rho(i,j)*((cv*T(i,j))+0.5*((u(i,j)^2)+(v(i,j)^2)))+p(i,j))*v(i,j)-u(i,j)*tau(i,j,1)-v(i,j)*tau(i,j,2)+q(i,j);
    end
end

end
