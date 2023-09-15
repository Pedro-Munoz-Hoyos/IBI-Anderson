% ==================================
% SUBSONIC-SUPERSONIC FLOW (CONSERVATIVE FORM)
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
U1=zeros(N,1);
U2=0.59*ones(N,1);
U3=zeros(N,1);

for i=1:N
    if (x(i)<=0.5)
        U1(i)=A(i);
        U3(i)=U1(i)*((1/(gam-1))+0.5*gam*((U2(i)/U1(i))^2));
    elseif (x(i)>0.5 && x(i)<=1.5)
        U1(i)=(1-0.366*(x(i)-0.5))*A(i);
        U3(i)=U1(i)*(((1-0.167*(x(i)-0.5))/(gam-1))+0.5*gam*((U2(i)/U1(i))^2));
    else
        U1(i)=(0.634-0.3879*(x(i)-1.5))*A(i);
        U3(i)=U1(i)*(((0.833-0.3507*(x(i)-1.5))/(gam-1))+0.5*gam*((U2(i)/U1(i))^2));
    end
end

% Residuals
r_U1=[]; 
r_U2=[];
r_U3=[];
nv=[];

n=1;
while(n<1401)
    % Time step
    dt=1000;
    for i=1:N
        dt_test=C*dx/(sqrt((gam-1)*((U3(i)/U1(i))-0.5*gam*((U2(i)/U1(i))^2)))+(U2(i)/U1(i)));
        if(dt_test<dt)
            dt=dt_test;
        end
    end

    r_U1(n)=0;
    r_U2(n)=0;
    r_U3(n)=0;
    for i=1:N-1
        % Predictor step  
        F1=U2(i);
        F11=U2(i+1);
        
        F2=((U2(i)^2)/U1(i))+((gam-1)/gam)*(U3(i)-0.5*gam*((U2(i)^2)/U1(i)));
        F21=((U2(i+1)^2)/U1(i+1))+((gam-1)/gam)*(U3(i+1)-0.5*gam*((U2(i+1)^2)/U1(i+1)));

        F3=gam*U2(i)*(U3(i)-(0.5*(gam-1)*((U2(i)^2)/U1(i))))/U1(i);
        F31=gam*U2(i+1)*(U3(i+1)-(0.5*(gam-1)*((U2(i+1)^2)/U1(i+1))))/U1(i+1);

        J2=(gam-1)*(U3(i)-0.5*gam*((U2(i)^2)/U1(i)))*(A(i+1)-A(i))/(A(i)*gam*dx);
        
        dU1dt_p=-(F11-F1)/dx;
        dU2dt_p=-((F21-F2)/dx)+J2;
        dU3dt_p=-(F31-F3)/dx;

        U1_m(i)=U1(i)+dU1dt_p*dt;
        U2_m(i)=U2(i)+dU2dt_p*dt;
        U3_m(i)=U3(i)+dU3dt_p*dt;

        if i>1
            % Corrector step
            F1=U2_m(i-1);
            F11=U2_m(i);
            
            F2=((U2_m(i-1)^2)/U1_m(i-1))+((gam-1)/gam)*(U3_m(i-1)-0.5*gam*((U2_m(i-1)^2)/U1_m(i-1)));
            F21=((U2_m(i)^2)/U1_m(i))+((gam-1)/gam)*(U3_m(i)-0.5*gam*((U2_m(i)^2)/U1_m(i)));
    
            F3=gam*U2_m(i-1)*(U3_m(i-1)-(0.5*(gam-1)*((U2_m(i-1)^2)/U1_m(i-1))))/U1_m(i-1);
            F31=gam*U2_m(i)*(U3_m(i)-(0.5*(gam-1)*((U2_m(i)^2)/U1_m(i))))/U1_m(i);
    
            J2=(gam-1)*(U3_m(i)-0.5*gam*((U2_m(i)^2)/U1_m(i)))*(A(i)-A(i-1))/(A(i)*gam*dx);

            dU1dt_c=-(F11-F1)/dx;
            dU2dt_c=-((F21-F2)/dx)+J2;
            dU3dt_c=-(F31-F3)/dx;

            dU1dt_av=(dU1dt_p+dU1dt_c)/2;
            dU2dt_av=(dU2dt_p+dU2dt_c)/2;
            dU3dt_av=(dU3dt_p+dU3dt_c)/2;

            if(dU1dt_av>r_U1(n))
                r_U1(n)=dU1dt_av;
            end
            if(dU2dt_av>r_U2(n))
                r_U2(n)=dU2dt_av;
            end
            if(dU3dt_av>r_U3(n))
                r_U3(n)=dU3dt_av;
            end

            U1_n1(i)=U1(i)+dU1dt_av*dt;
            U2_n1(i)=U2(i)+dU2dt_av*dt;
            U3_n1(i)=U3(i)+dU3dt_av*dt;
        end

    end
   
    % Subsonic inflow B.Cs.
    U1_n1(1)=A(1);
    U2_n1(1)=2*U2_n1(2)-U2_n1(3);
    U3_n1(1)=U1_n1(1)*((1/(gam-1))+0.5*gam*((U2_n1(1)/U1_n1(1))^2));
    
    % Supersonic outflow B.Cs.
    U1_n1(N)=2*U1_n1(N-1)-U1_n1(N-2);
    U2_n1(N)=2*U2_n1(N-1)-U2_n1(N-2);
    U3_n1(N)=2*U3_n1(N-1)-U3_n1(N-2);

    U1=U1_n1;
    U2=U2_n1;
    U3=U3_n1;
    
    nv(n)=n;
    fig1=figure(1);
    refresh(fig1);
    set(gcf,'color','white');
    plot(nv,r_U1);
    hold on;
    plot(nv,r_U2);
    hold on;
    plot(nv,r_U3);
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
rho_n=zeros(N,1);
V_n=zeros(N,1);
T_n=zeros(N,1);
p_n=zeros(N,1);
M_n=zeros(N,1);
mdot_n=zeros(N,1);
for i=1:N
    rho_n(i)=U1(i)/A(i);
    V_n(i)=U2(i)/U1(i);
    T_n(i)=(gam-1)*((U3(i)/U1(i))-0.5*gam*((U2(i)/U1(i))^2));
    p_n(i)=rho_n(i)*T_n(i);
    M_n(i)=V_n(i)/sqrt(T_n(i));
    mdot_n(i)=rho_n(i)*V_n(i)*A(i);
end

Comparator=zeros(N,11);
Comparator(:,1)=round(x,3);
Comparator(:,2)=round(A,3);
Comparator(:,3)=round(rho_n,3);
Comparator(:,4)=round(V_n,3);
Comparator(:,5)=round(T_n,3);
Comparator(:,6)=round(p_n,3);
Comparator(:,7)=round(M_n,3);
Comparator(:,8)=round(mdot_n,3);
Comparator(:,9)=round(U1,3);
Comparator(:,10)=round(U2,3);
Comparator(:,11)=round(U3,3);
