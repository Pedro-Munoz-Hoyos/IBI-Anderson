% ==================================
% COUETTE FLOW
% ==================================
close all;
% ==================================
% SIMULATION PARAMETERS
% ==================================
ReD=5000;   % Reynolds number based on the plate's distance
N=21;       % Number of y-grid points
conv=10^-6; % Convergence criterion

%% ==================================
% NUMERICAL ENGINE
% ==================================
% Initial conditions
t=0;
dy=1/(N-1);
dt=0.5*ReD*(dy^2);
y=0:dy:1;
u=zeros(N,1);
u(N)=1;

c=0;
tv=[];
resv=[];
res=conv+1;
dudtp=zeros(N-2,1);
um=zeros(N,1);
while(res>conv)   
    % Predictor step
    for j=1:N
        if j==1
            dudtp(j)=(2*u(j)-5*u(j+1)+4*u(j+2)-u(j+3))/(ReD*(dy^2));
        elseif j==N
            dudtp(j)=(2*u(j)-5*u(j-1)+4*u(j-2)-u(j-3))/(ReD*(dy^2));
        else
            dudtp(j)=(u(j-1)-2*u(j)+u(j+1))/(ReD*(dy^2));
        end
        um(j)=u(j)+dudtp(j)*dt;
    end
    
    res=0;
    % Corrector step
    for j=2:N-1
        dudtc=(um(j-1)-2*um(j)+um(j+1))/(ReD*(dy^2));
        rest=0.5*(dudtc+dudtp(j))*dt;
        u(j)=u(j)+rest;
        
        % Residual's update
        if abs(rest)>res
            res=abs(rest);
        end
    end

    c=c+1;
    t=t+dt;
    tv(c)=t;
    resv(c)=res;

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
