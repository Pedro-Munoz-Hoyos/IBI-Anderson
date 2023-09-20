function [pn1,un1,vn1] = simple(us,vs,ps,rho,mu,dx,dy,dt)
% Function in charge of implementing the pressure correction method by
% means of the SIMPLE (Semi-Implicit Pressure Linked-equations) algorithm.
% The imposed boundary conditions are:
% - Inflow:  Dirichlet (p=p1, v=0)
% - Walls:   Dirichlet (u=0, v=0)
%            Neumann   (dp/dy=0, dp'/dy=0)
% - Outflow: Dirichlet (p=p2)

nux=size(us,1);
nuy=size(us,2);
nvx=size(vs,1);
nvy=size(vs,2);
npx=size(ps,1);
npy=size(ps,2);

usn1=zeros(nux,nuy);
vsn1=zeros(nvx,nvy);
un1=zeros(nux,nuy);
vn1=zeros(nvx,nvy);
pn1=zeros(npx,npy);

% x-momentum equation
for i=2:(nux-1)
    for j=2:(nuy-1)
        vm=0.5*(vs(i,j+1)+vs(i-1,j+1));
        vmm=0.5*(vs(i,j)+vs(i-1,j));
        A1=-0.5*rho*((((us(i+1,j)^2)-(us(i-1,j)^2))/dx)+(((us(i,j+1)*vm)-(us(i,j-1)*vmm))/dy));
        A2=mu*(((us(i+1,j)-2*us(i,j)+us(i-1,j))/(dx^2))+((us(i,j+1)-2*us(i,j)+us(i,j-1))/(dy^2)));    
        usn1(i,j)=us(i,j)+((A1+A2)*dt-(dt*(ps(i,j)-ps(i-1,j))/dx))/rho;
    end
end

% y-momentum equation
for i=2:(nvx-1)
    for j=2:(nvy-1)
        um=0.5*(us(i+1,j-1)+us(i+1,j));
        umm=0.5*(us(i,j-1)+us(i,j));
        B1=-0.5*rho*((((vs(i,j+1)^2)-(vs(i,j-1)^2))/dy)+((vs(i+1,j)*um-vs(i-1,j)*umm)/dx));
        B2=mu*(((vs(i+1,j)-2*vs(i,j)+vs(i-1,j))/(dx^2))+((vs(i,j+1)-2*vs(i,j)+vs(i,j-1))/(dy^2)));
        vsn1(i,j)=vs(i,j)+((B1+B2)*dt-(dt*(ps(i,j)-ps(i,j-1))/dy))/rho;
    end
end

% Relaxation technique with overrelaxation factor (w)
w=1.5;
err=1;
b=-dt/(dx^2);
c=-dt/(dy^2);
a=-2*(b+c);
pp=zeros(npx,npy);
while(err>10^-10)
    err=0;
    for i=2:npx-1
        for j=2:npy-1
            prev=pp(i,j);
            d=rho*(((usn1(i+1,j)-usn1(i,j))/dx)+((vsn1(i,j+1)-vsn1(i,j))/dy));
            ppmean=-(b*(pp(i+1,j)+pp(i-1,j))+c*(pp(i,j+1)+pp(i,j-1))+d)/a;
            pp(i,j)=prev+w*(ppmean-prev);

            % Relaxation technique's boundary conditions (dp'/dy)_w=0 and
            % p'_1=p'_Nx=0
            if(j==2)
                pp(i,1)=pp(i,2);
            elseif(j==npy-1)
                pp(i,npy)=pp(i,npy-1);
            end

            % Check convergence of the relaxation technique
            errt=abs(pp(i,j)-prev);
            if(errt>err)
                err=errt;
            end
        end
    end
end

% Corrected pressure (equivalent to the next iteration's guessed value)
ap=0.1; % Underrelaxation factor
for i=1:npx
    for j=1:npy
        pn1(i,j)=ps(i,j)+ap*pp(i,j);
    end
end

% Corrected x-velocity (equivalent to the next iteration's guessed value)
for i=2:(nux-1)
    for j=2:(nuy-1)
        un1(i,j)=usn1(i,j)-(dt*(pp(i,j)-pp(i-1,j))/(dx*rho));
    end
end

% Corrected y-velocity (equivalent to the next iteration's guessed value)
for i=2:(nvx-1)
    for j=2:(nvy-1)
        vn1(i,j)=vsn1(i,j)-(dt*(pp(i,j)-pp(i,j-1))/(dy*rho));
    end
end

% u-boundary conditions
un1(1,:)=un1(2,:);          % Inlet
un1(:,1)=us(:,1);           % Lower wall     
un1(:,nuy)=us(:,nuy);       % Upper wall
un1(nux,:)=un1(nux-1,:);    % Outlet

% v-boundary conditions
vn1(1,:)=0;                 % Inlet
vn1(:,1)=-vn1(:,2);         % Lower wall     
vn1(:,nvy)=-vn1(:,nvy-1);   % Upper wall
vn1(nvx,:)=vn1(nvx-1,:);    % Outlet

end
