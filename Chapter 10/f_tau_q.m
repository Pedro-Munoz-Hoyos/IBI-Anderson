function [tau,q] = tau_q(mu,cp,u,v,T,dx,dy,pred,E)
% Function in charge of obtaining the shear and normal stresses together
% with the thermal conduction term, taking into account whether it is being
% used in the predictor or corrector step in order to mantain second-order
% accuracy.

nx=size(u,1);
ny=size(u,2);
tau=zeros(nx,ny,2);
q=zeros(nx,ny);

if (pred==1)
    for i=1:nx
        for j=1:ny  
            if(E==1)
                if i==1
                    dudx=(u(i+1,j)-u(i,j))/dx;
                    dvdx=(v(i+1,j)-v(i,j))/dx;
                    dTdx=(T(i+1,j)-T(i,j))/dx;
                else
                    dudx=(u(i,j)-u(i-1,j))/dx;
                    dvdx=(v(i,j)-v(i-1,j))/dx;
                    dTdx=(T(i,j)-T(i-1,j))/dx;
                end
                if j==1
                    dudy=(u(i,j+1)-u(i,j))/dy;
                    dvdy=(v(i,j+1)-v(i,j))/dy;
                elseif j==ny
                    dudy=(u(i,j)-u(i,j-1))/dy;
                    dvdy=(v(i,j)-v(i,j-1))/dy;
                else
                    dudy=0.5*(u(i,j+1)-u(i,j-1))/dy;
                    dvdy=0.5*(v(i,j+1)-v(i,j-1))/dy;
                end
                div=dudx+dvdy;
                tau(i,j,1)=mu(i,j)*(dudy+dvdx);
                tau(i,j,2)=2*mu(i,j)*((-div/3)+dudx);
                k=mu(i,j)*cp/0.71;
                q(i,j)=-k*dTdx;
            else
                if j==1
                    dudy=(u(i,j+1)-u(i,j))/dy;
                    dvdy=(v(i,j+1)-v(i,j))/dy;
                    dTdy=(T(i,j+1)-T(i,j))/dy;
                else
                    dudy=(u(i,j)-u(i,j-1))/dy;
                    dvdy=(v(i,j)-v(i,j-1))/dy;
                    dTdy=(T(i,j)-T(i,j-1))/dy;
                end
                if i==1
                    dudx=(u(i+1,j)-u(i,j))/dx;
                    dvdx=(v(i+1,j)-v(i,j))/dx;
                elseif i==nx
                    dudx=(u(i,j)-u(i-1,j))/dx;
                    dvdx=(v(i,j)-v(i-1,j))/dx;
                else
                    dudx=0.5*(u(i+1,j)-u(i-1,j))/dx;
                    dvdx=0.5*(v(i+1,j)-v(i-1,j))/dx;
                end
                div=dudx+dvdy;
                tau(i,j,1)=mu(i,j)*(dudy+dvdx);
                tau(i,j,2)=2*mu(i,j)*((-div/3)+dvdy);
                k=mu(i,j)*cp/0.71;
                q(i,j)=-k*dTdy;
            end
        end
    end
else
    for i=1:nx
        for j=1:ny  
            if(E==1)
                if i==nx
                    dudx=(u(i,j)-u(i-1,j))/dx;
                    dvdx=(v(i,j)-v(i-1,j))/dx;
                    dTdx=(T(i,j)-T(i-1,j))/dx;
                    
                else
                    dudx=(u(i+1,j)-u(i,j))/dx;
                    dvdx=(v(i+1,j)-v(i,j))/dx;
                    dTdx=(T(i+1,j)-T(i,j))/dx;
                end
                if j==1
                    dudy=(u(i,j+1)-u(i,j))/dy;
                    dvdy=(v(i,j+1)-v(i,j))/dy;
                elseif j==ny
                    dudy=(u(i,j)-u(i,j-1))/dy;
                    dvdy=(v(i,j)-v(i,j-1))/dy;
                else
                    dudy=0.5*(u(i,j+1)-u(i,j-1))/dy;
                    dvdy=0.5*(v(i,j+1)-v(i,j-1))/dy;
                end
                div=dudx+dvdy;
                tau(i,j,1)=mu(i,j)*(dudy+dvdx);
                tau(i,j,2)=2*mu(i,j)*((-div/3)+dudx);
                k=mu(i,j)*cp/0.71;
                q(i,j)=-k*dTdx;
            else
                if j==ny
                    dudy=(u(i,j)-u(i,j-1))/dy;
                    dvdy=(v(i,j)-v(i,j-1))/dy;
                    dTdy=(T(i,j)-T(i,j-1))/dy;
                else
                    dudy=(u(i,j+1)-u(i,j))/dy;
                    dvdy=(v(i,j+1)-v(i,j))/dy;
                    dTdy=(T(i,j+1)-T(i,j))/dy;
                end
                if i==1
                    dudx=(u(i+1,j)-u(i,j))/dx;
                    dvdx=(v(i+1,j)-v(i,j))/dx;
                elseif i==nx
                    dudx=(u(i,j)-u(i-1,j))/dx;
                    dvdx=(v(i,j)-v(i-1,j))/dx;
                else
                    dudx=0.5*(u(i+1,j)-u(i-1,j))/dx;
                    dvdx=0.5*(v(i+1,j)-v(i-1,j))/dx;
                end
                div=dudx+dvdy;
                tau(i,j,1)=mu(i,j)*(dudy+dvdx);
                tau(i,j,2)=2*mu(i,j)*((-div/3)+dvdy);
                k=mu(i,j)*cp/0.71;
                q(i,j)=-k*dTdy;
            end
        end
    end
end

end
