% ==================================
% COUETTE FLOW
% ==================================
close all;
% ==================================
% SIMULATION PARAMETERS
% ==================================
ReD=5000;   % Reynolds number based on the plate's distance
E=1;                % Marching step related parameter
N=21;               % Number of y-grid points
conv=10^-8;         % Convergence criterion


%% ==================================
% NUMERICAL ENGINE
% ==================================
% Initial conditions
t=0;
dy=1/(N-1);
dt=E*ReD*(dy^2);
y=0:dy:1;
u=zeros(N,1);
u(N)=1;

c=1;
res=conv+1;
tv=[];
resv=[];
while(res>conv)
    [u,res]=talgorithm(u,E,dt);
    tv(c)=t;
    resv(c)=res;
    
    c=c+1;
    t=t+dt;
    fig1=figure(1);
    refresh(fig1);
    set(gcf,'color','white');
    plot(tv,resv);
    xlabel("t'");
    ylabel("max(\partial/\partialt')");
    title("Steady-state convergence");

    fig2=figure(2);
    refresh(fig2);
    set(gcf,'color','white');
    plot(u,y);
    xlabel("y'");
    ylabel("U'(y')");
    title("Nondimensional velocity distribution");
end
