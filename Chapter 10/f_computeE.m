function [E] = computeE(rho,mu,cv,cp,dx,dy,u,v,T,p,pred)
% Function in charge of obtaining the flux vector E from the primitive
% variables, taking into account whether it is being calculated in the 
% predictor or corrector step in order to mantain second-order accuracy.


nx=size(u,1);
ny=size(u,2);
E=zeros(nx,ny,4);
[tau,q]=tau_q(mu,cp,u,v,T,dx,dy,pred,1);

for i=1:nx
    for j=1:ny
        E(i,j,1)=rho(i,j)*u(i,j);
        E(i,j,2)=E(i,j,1)*u(i,j)+p(i,j)-tau(i,j,2);
        E(i,j,3)=E(i,j,1)*v(i,j)-tau(i,j,1);
        E(i,j,4)=(rho(i,j)*((cv*T(i,j))+0.5*((u(i,j)^2)+(v(i,j)^2)))+p(i,j))*u(i,j)-u(i,j)*tau(i,j,2)-v(i,j)*tau(i,j,1)+q(i,j);
    end
end

end
