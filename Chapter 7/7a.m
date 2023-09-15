% ==================================
% SUBSONIC-SUPERSONIC FLOW
% ==================================
% FLOW PARAMETERS
% ==================================
N=31;    % Number of axial grid points (Odd number)
C=0.5;   % Courant Number
gam=1.4; % Specific heats' ratio


% ==================================
% GRID GENERATION
% ==================================
A=zeros(N,1); % Non-dimensional Area
dx=3/(N-1);   % Non-dimensional x-step
x=0:dx:3;
for i=1:N
    A(i)=1+2.2*((x(i)-1.5)^2);
end        


% ==================================
% NUMERICAL ENGINE
% ==================================
% Initial conditions
rho_n=1-0.3146.*x;
T_n=1-0.2314.*x;
V_n=(0.1+1.09.*x).*(T_n.^0.5);

% Residuals
r_rho=[]; 
r_V=[];
r_T=[];
nv=[];

n=1;
while(n<1401)
    % Time step
    dt=1000;
    for i=1:N
        dt_test=C*dx/(sqrt(T_n(i))+V_n(i));
        if(dt_test<dt)
            dt=dt_test;
        end
    end

    r_rho(n)=0;
    r_V(n)=0;
    r_T(n)=0;
    for i=1:N-1
        % Predictor step
        dA=log(A(i+1))-log(A(i));
        drho=rho_n(i+1)-rho_n(i);
        dV=V_n(i+1)-V_n(i);
        dT=T_n(i+1)-T_n(i);
        
        dVdt_p=(-V_n(i)*dV-(1/gam)*(dT+(T_n(i)/rho_n(i))*drho))/dx;
        if i==1
            drhodt_p=0;
            dTdt_p=0;
        else
            drhodt_p=(-rho_n(i)*dV-rho_n(i)*V_n(i)*dA-V_n(i)*drho)/dx;
            dTdt_p=(-V_n(i)*dT-(gam-1)*T_n(i)*(dV+V_n(i)*dA))/dx;
        end
        
        rho_m(i)=rho_n(i)+drhodt_p*dt;
        V_m(i)=V_n(i)+dVdt_p*dt;
        T_m(i)=T_n(i)+dTdt_p*dt;

        if i>1
            % Corrector step
            dA=log(A(i))-log(A(i-1));
            drho=rho_m(i)-rho_m(i-1);
            dV=V_m(i)-V_m(i-1);
            dT=T_m(i)-T_m(i-1);
            
            drhodt_c=(-rho_m(i)*dV-rho_m(i)*V_m(i)*dA-V_m(i)*drho)/dx;
            dVdt_c=(-V_m(i)*dV-(1/gam)*(dT+(T_m(i)/rho_m(i))*drho))/dx;
            dTdt_c=(-V_m(i)*dT-(gam-1)*T_m(i)*(dV+V_m(i)*dA))/dx;
            
            drhodt_av=(drhodt_p+drhodt_c)/2;
            dVdt_av=(dVdt_p+dVdt_c)/2;
            dTdt_av=(dTdt_p+dTdt_c)/2;

            if(drhodt_av>r_rho(n))
                r_rho(n)=drhodt_av;
            end
            if(dVdt_av>r_V(n))
                r_V(n)=dVdt_av;
            end
            if(dTdt_av>r_T(n))
                r_T(n)=dTdt_av;
            end

            rho_n1(i)=rho_n(i)+drhodt_av*dt;
            V_n1(i)=V_n(i)+dVdt_av*dt;
            T_n1(i)=T_n(i)+dTdt_av*dt;
        end

    end
   
    % Subsonic inflow B.Cs.
    rho_n1(1)=1;
    V_n1(1)=2*V_n1(2)-V_n1(3);
    T_n1(1)=1;
    
    % Supersonic outflow B.Cs.
    rho_n1(N)=2*rho_n1(N-1)-rho_n1(N-2);
    V_n1(N)=2*V_n1(N-1)-V_n1(N-2);
    T_n1(N)=2*T_n1(N-1)-T_n1(N-2);

    rho_n=rho_n1;
    V_n=V_n1;
    T_n=T_n1;
    
    nv(n)=n;
    fig1=figure(1);
    refresh(fig1);
    set(gcf,'color','white');
    plot(nv,r_rho);
    hold on;
    plot(nv,r_V);
    hold on;
    plot(nv,r_T);
    hold off;
    xlabel("Iterations");
    ylabel("max(\partial/\partialt)");
    title("Steady-state convergence");
    pause(1/100);

    n=n+1;
end


% ==================================
% POST-PROCESSING
% ==================================
p_n=zeros(N,1);
M_n=zeros(N,1);
mdot_n=zeros(N,1);
for i=1:N
    p_n(i)=rho_n(i)*T_n(i);
    M_n(i)=V_n(i)/sqrt(T_n(i));
    mdot_n(i)=rho_n(i)*V_n(i)*A(i);
end

Comparator=zeros(N,8);
Comparator(:,1)=round(x,3);
Comparator(:,2)=round(A,3);
Comparator(:,3)=round(rho_n,3);
Comparator(:,4)=round(V_n,3);
Comparator(:,5)=round(T_n,3);
Comparator(:,6)=round(p_n,3);
Comparator(:,7)=round(M_n,3);
Comparator(:,8)=round(mdot_n,3);
