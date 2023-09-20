% ==================================
% COUETTE FLOW (PRESSURE CORRECTION METHOD)
% ==================================
close all;
% ==================================
% FLOW PARAMETERS
% ==================================
uu=0;            % Upper plate's x-velocity [ft/s]
ud=0;            % Bottom plate's x-velocity [ft/s]
p1=14.6959;      % Inlet's pressure [psi] 
p2=0.99*p1;      % Outelt's pressure [psi] 
rho=0.002377;    % Air density [slug/ft3]
mu=3.737*10^(-7);% Air kinematic viscosity [slug/(ft·s)]


%% ==================================
% SIMULATION PARAMETERS
% ==================================
L=0.5;           % Length of the domain [ft]
H=0.01;          % Height of the domain [ft]
Nx=21;           % Number of x-grid points
Ny=11;           % Number of y-grid points
dt=0.001;        % Time step [s]


%% ==================================
% NUMERICAL ENGINE
% ==================================
p=zeros(Nx,Ny);
u=zeros(Nx+1,Ny);
v=zeros(Nx+1,Ny+1);
dx=L/(Nx-1);
dy=H/(Ny-1);

% Initial conditions
p(1,:)=p1;
p(Nx,:)=p2;
u(:,1)=ud;
u(:,Ny)=uu;

c=0;
cv=[];
res=[];
while(c<10000)
    [p,u,v]=simple(u,v,p,rho,mu,dx,dy,dt);
    c=c+1;
    fprintf("Iteration: %i\n",c);
end


%% ==================================
% POST-PROCESSING
% ==================================
x=zeros(Nx,Ny);
y=zeros(Nx,Ny);
uij=zeros(Nx,Ny);
vij=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        x(i,j)=(i-1)*dx;
        y(i,j)=(j-1)*dy;
        uij(i,j)=0.5*(u(i+1,j)+u(i,j));
        vij(i,j)=0.5*(v(i,j+1)+v(i,j));
    end
end
contourf(x,y,uij,80,'LineColor','none');
set(gcf,'color','white');
xlabel("x[m]");
ylabel("y[m]");
title("X-velocity distribution U(x,y)");
