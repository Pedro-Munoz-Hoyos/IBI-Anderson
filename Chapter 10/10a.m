% ==================================
% SUPERSONIC FLOW OVER A FLAT PLATE (N-S)
% ==================================
close all;
% - A FIXED WALL TEMPERATURE B.C. IS CURRENTLY IMPLEMENTED
% - CHANGE THE "bcs" FUNCTION TO IMPLEMENT AN ADIABATIC WALL B.C.
%   (FOR AN ADIABATIC WALL CASE, USE K=0.5)
% ==================================
% FLOW PARAMETERS
% ==================================
gam=1.4;               % Specific heats' ratio
R=287.05;              % Specific air gas constant [J/(Kg·K)]
cv=R/(gam-1);          % Specific heat at constant volume [J/(Kg·K)]
cp=gam*cv;             % Specific heat at constant pressure [J/(Kg·K)]
M_fs=4;                % Free-stream Mach Number
p_fs=101325;           % Free-stream pressure [Pa]
T_fs=288.16;           % Free-stream temperature [K]
rho_fs=p_fs/(R*T_fs);  % Free-stream density [kg/m3]
c_fs=sqrt(gam*R*T_fs); % Free-stream speed of sound [m/s]
mu_fs=1.7894*10^(-5);  % Free-stream kinematic viscosity [kg/(m·s)]
T_w=T_fs;              % Wall temperature [K]



%% ==================================
% SIMULATION PARAMETERS
% ==================================
L=0.00001;                             % Length of the domain [m]
H=25*L/sqrt(rho_fs*M_fs*c_fs*L/mu_fs); % Height of the domain [m]
Nx=70;                                 % Number of x-grid points
Ny=70;                                 % Number of y-grid points
K=0.6;                                 % Courant number


%% ==================================
% NUMERICAL ENGINE
% ==================================
dx=L/(Nx-1);
dy=H/(Ny-1);

% Initial conditions
rho=rho_fs*ones(Nx,Ny);
u=M_fs*c_fs*ones(Nx,Ny);
v=zeros(Nx,Ny);
T=T_fs*ones(Nx,Ny);
p=p_fs*ones(Nx,Ny);
u(:,1)=zeros(Nx,1);
T(2:Nx,1)=T_w*ones(Nx-1,1);

mu=computemu(mu_fs,T_fs,T);
U=computeU(rho,cv,u,v,T);
E=computeE(rho,mu,cv,cp,dx,dy,u,v,T,p,1);
F=computeF(rho,mu,cv,cp,dx,dy,u,v,T,p,1);

c=0;
Um=U;
dudtp=zeros(Nx-2,Ny-2,4);
while(c<9000)
    % TIME STEP
    dt=computedt(rho,mu,u,v,T,gam,R,dx,dy,K);

    % PREDICTOR STEP
    for i=2:Nx-1
        for j=2:Ny-1
            for b=1:4
                dudtp(i-1,j-1,b)=-((E(i+1,j,b)-E(i,j,b))/dx)-((F(i,j+1,b)-F(i,j,b))/dy);
                Um(i,j,b)=U(i,j,b)+dudtp(i-1,j-1,b)*dt;
            end
        end
    end  
    [rhom,um,vm,Tm,pm]=decodevars(Um,cv,R);
    [rhom,um,vm,Tm,pm,Um(2:Nx,1,:),Um(Nx,2:Ny-1,:)]=bcs(um,vm,Tm,pm,rhom,cv,R);

    % CORRECTOR STEP
    mum=computemu(mu_fs,T_fs,Tm);
    Em=computeE(rhom,mum,cv,cp,dx,dy,um,vm,Tm,pm,0);
    Fm=computeF(rhom,mum,cv,cp,dx,dy,um,vm,Tm,pm,0);
    for i=2:Nx-1
        for j=2:Ny-1
            for b=1:4
                dudtc=-((Em(i,j,b)-Em(i-1,j,b))/dx)-((Fm(i,j,b)-Fm(i,j-1,b))/dy);
                U(i,j,b)=U(i,j,b)+0.5*(dudtc+dudtp(i-1,j-1,b))*dt;
            end
        end
    end
    [rho,u,v,T,p]=decodevars(U,cv,R);
    [rho,u,v,T,p,U(2:Nx,1,:),U(Nx,2:Ny-1,:)]=bcs(u,v,T,p,rho,cv,R);

    mu=computemu(mu_fs,T_fs,T);
    E=computeE(rho,mu,cv,cp,dx,dy,u,v,T,p,1);
    F=computeF(rho,mu,cv,cp,dx,dy,u,v,T,p,1);

    c=c+1;
    fprintf("c:%i\n",c);
end


%% ==================================
% POST-PROCESSING
% ==================================
x=0:dx:L;
y=0:dy:H;
ynorm=zeros(Ny,1);
M=zeros(Ny,1);
for j=1:Ny
    ynorm(j)=y(j)*sqrt(rho(Nx,j)*u(Nx,j)*L/mu(Nx,j))/L;
    M(j)=sqrt(((u(Nx,j)^2)+(v(Nx,j)^2))/(gam*R*T(Nx,j)));
end

figure(1);
set(gcf,'color','white');
plot(x,p(:,1)/p_fs);
xlabel("x[m]");
ylabel("p/p_{\infty}");
title("Normalized surface pressure distributions");


figure(2);
set(gcf,'color','white');
plot(p(Nx,:)/p_fs,ynorm);
xlabel("p/p_{\infty}");
ylabel("y_{norm}");
title("Normalized pressure profiles (at t.e.)");

figure(3);
set(gcf,'color','white');
plot(T(Nx,:)/T_fs,ynorm);
xlabel("T/T_{\infty}");
ylabel("y_{norm}");
title("Normalized temperature profiles (at t.e.)");

figure(4);
set(gcf,'color','white');
plot(u(Nx,:)/(M_fs*c_fs),ynorm);
xlabel("u/u_{\infty}");
ylabel("y_{norm}");
title("Normalized x-velocity profiles (at t.e.)");

figure(5);
set(gcf,'color','white');
plot(M,ynorm);
xlabel("M/M_{\infty}");
ylabel("y_{norm}");
title("Normalized Mach number profiles (at t.e.)");

figure(6);
set(gcf,'color','white');
contourf(transpose(p),80,'LineColor','None');
xlabel("x[m]");
ylabel("y[m]");
title("Pressure distribution p(x,y)");
