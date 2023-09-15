% ==================================
% PRANDTL-MEYER EXPANSION WAVE
% ==================================
close all;
% ==================================
% FLOW PARAMETERS
% ==================================
M1=2;               % Inflow Mach number
p1=1.01*10^5;       % Inflow pressure [Pa]
rho1=1.23;          % Inflow density  [kg/m3]
T1=286.1;           % Inflow temperature [K]
gam=1.4;            % Specific heats' ratio
R=287.052874247;    % Specific air gas constant [J/(Kg·K)]


%% ==================================
% SIMULATION PARAMETERS
% ==================================
L=67;               % Length of the domain [m]
E=10;               % X-location of the expansion corner [m]
H=40;               % Upper boundary height [m]
tht=5.352*pi/180;   % Expansion angle [rad]
C=0.5;              % Courant Number
Cy=0.6;             % Artificial viscosity constant
Ny=81;              % Number of y-grid points


%% ==================================
% GRID GENERATION
% ==================================
% Computational space
xi=zeros(2,1);
deta=1/(Ny-1);
eta=0:deta:1;

% Physical space
x=zeros(Ny,1);
y=zeros(Ny,1);
% Note that the matter of grid generation for the case of a downstream
% marching solution depends upon the actual size of the marching step. As a
% result, the grid is generated during the actual numerical solution.


%% ==================================
% NUMERICAL ENGINE
% ==================================
% Initial conditions
rho=rho1*ones(Ny,1);
u=M1*sqrt(gam*R*T1)*ones(Ny,1);
v=zeros(Ny,1);
p=p1*ones(Ny,1);
T=T1*ones(Ny,1);
M=M1*ones(Ny,1);
F(:,1)=rho1*u;
F(:,2)=F(:,1)*u(1)+p1*ones(Ny,1);
F(:,3)=zeros(Ny,1);
F(:,4)=(gam*p1*u(1)/(gam-1))*ones(Ny,1)+F(:,1)*0.5*(u(1)^2);
G(:,1)=zeros(Ny,1);
G(:,2)=zeros(Ny,1);
G(:,3)=p1*ones(Ny,1);
G(:,4)=zeros(Ny,1);

n=1;
dFp=zeros(Ny,4);
Fm=zeros(Ny,4);
rhom=zeros(Ny,1);
Gm=zeros(Ny,4);
while (xi(n)<L)
    % Corner criterion
    if(xi(n)<=E)
        h=H;
        A=0;
        wang=0;
    else
        h=(H+(xi(n)-E)*tan(tht));
        A=tan(tht)/h;
        wang=tht;
    end

    % Physical space generation
    for j=1:Ny
        x(Ny-j+1,n)=xi(n);
        y(Ny-j+1,n)=H-h*deta*(j-1);
    end

    % Marching step
    [dxi]=step(C,h,u,v,M);

    % Predictor step
    for j=1:Ny      
        % Artificial viscosity
        if(j>1 && j<Ny)
            [S]=artivis(F(j-1,:),F(j,:),F(j+1,:),rho(j,n),Cy);
        else
            S=zeros(4,1);
        end
        % Predicted values
        for i=1:4
            if(j==Ny)
                % Second order backward difference
                dFp(j,i)=((A*(1-eta(j))*(4*F(j-1,i)-3*F(j,i)-F(j-2,i)))+((4*G(j-1,i)-3*G(j,i)-G(j-2,i))/h))/(2*deta);
            else
                dFp(j,i)=((A*(1-eta(j))*(F(j,i)-F(j+1,i)))+((G(j,i)-G(j+1,i))/h))/deta;
            end
            Fm(j,i)=F(j,i)+dFp(j,i)*dxi+S(i);
        end
        [rhom(j)]=decoderho(Fm(j,:),gam);
        [Gm(j,:)]=decodeg(Fm(j,:),rhom(j),gam);
    end

    % Corrector step
    for j=1:Ny
        % Artificial viscosity
        if(j>1 && j<Ny)
            [S]=artivis(Fm(j-1,:),Fm(j,:),Fm(j+1,:),rhom(j),Cy);
        else
            S=zeros(4,1);
        end
        % Corrected values
        for i=1:4
            if(j==1)
                % Second order forward difference
                dFc=((A*(1-eta(j))*(3*Fm(j,i)-4*Fm(j+1,i)+Fm(j+2,i)))+((3*Gm(j,i)-4*Gm(j+1,i)+Gm(j+2,i))/h))/(2*deta);
            else
                dFc=((A*(1-eta(j))*(Fm(j-1,i)-Fm(j,i)))+((Gm(j-1,i)-Gm(j,i))/h))/deta;
            end
            F(j,i)=F(j,i)+0.5*(dFc+dFp(j,i))*dxi+S(i);
        end
        [rho(j,n+1)]=decoderho(F(j,:),gam);
        [G(j,:)]=decodeg(F(j,:),rho(j,n+1),gam);
        [u(j,n+1),v(j,n+1),p(j,n+1),T(j,n+1),M(j,n+1)]=decodevars(F(j,:),rho(j,n+1),gam,R);
    end
    

    % Boundary conditions
    [rho(1,n+1),u(1,n+1),v(1,n+1),p(1,n+1),T(1,n+1),M(1,n+1),F(1,:),G(1,:)] = bcs(u(1,n+1),v(1,n+1),p(1,n+1),T(1,n+1),M(1,n+1),gam,R,wang);
    [rho(Ny,n+1),u(Ny,n+1),v(Ny,n+1),p(Ny,n+1),T(Ny,n+1),M(Ny,n+1),F(Ny,:),G(Ny,:)] = bcs(u(Ny,n+1),v(Ny,n+1),p(Ny,n+1),T(Ny,n+1),M(Ny,n+1),gam,R,0);

    xi(n+1)=xi(n)+dxi;
    fprintf("x:%f\n",xi(n+1));
    n=n+1;
end


%% ==================================
% POST-PROCESSING
% ==================================
Comparator=zeros(Ny,12);
Comparator(:,1)=round(y(:,n-1),3);
Comparator(:,2)=round(eta,3);
Comparator(:,3)=round(u(:,n-1),3);
Comparator(:,4)=round(v(:,n-1),3);
Comparator(:,5)=round(rho(:,n-1),3);
Comparator(:,6)=round(p(:,n-1),3);
Comparator(:,7)=round(T(:,n-1),3);
Comparator(:,8)=round(M(:,n-1),3);
Comparator(:,9)=round(F(:,1),3);
Comparator(:,10)=round(F(:,2),3);
Comparator(:,11)=round(F(:,3),3);
Comparator(:,12)=round(F(:,4),3);

% Trailing edge analytical solution
M2=M1;
err=1;
f1=sqrt((gam+1)/(gam-1))*atan(sqrt((gam-1)*((M1^2)-1)/(gam+1)))-atan(sqrt((M1^2)-1));
f2=f1+tht;
while err>10^-7
    err=abs(sqrt((gam+1)/(gam-1))*atan(sqrt((gam-1)*((M2^2)-1)/(gam+1)))-atan(sqrt((M2^2)-1))-f2);
    M2=M2+10^-7;
end

lwave=zeros(100,2);
twave=zeros(100,2);
for i=1:100
    lwave(i,1)=E+((xi(n-1)-E)/99)*(i-1);
    lwave(i,2)=((xi(n-1)-E)/99)*(i-1)*tan(asin(1/M1));
    twave(i,1)=E+((xi(n-1)-E)/99)*(i-1);
    twave(i,2)=((xi(n-1)-E)/99)*(i-1)*tan(asin(1/M2)-tht);
end
contourf(x,y,M(:,1:n-1),80,'LineColor','none');
hold on;
plot(lwave(:,1),lwave(:,2),'--k');
hold on;
plot(twave(:,1),twave(:,2),'--k');
hold off;
set(gcf,'color','white');
xlabel("x[m]");
ylabel("y[m]");
title("Mach number distribution M(x,y)");
