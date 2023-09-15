% ==================================
% SUBSONIC-SUPERSONIC FLOW (SHOCK CAPTURING)
% ==================================
% FLOW PARAMETERS
% ==================================
N=61;        % Number of axial grid points (Odd number)
C=0.5;       % Courant Number
gam=1.4;     % Specific heats' ratio
pepo=0.6784; % Pressure ratio across the nozzle
Cx=0.2;      % Artificial viscosity constant


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
    elseif (x(i)>1.5 && x(i)<=2.1)
        U1(i)=(0.634-0.702*(x(i)-1.5))*A(i);
        U3(i)=U1(i)*(((0.833-0.4908*(x(i)-1.5))/(gam-1))+0.5*gam*((U2(i)/U1(i))^2));
    else
        U1(i)=(0.5892+0.10228*(x(i)-2.1))*A(i);
        U3(i)=U1(i)*(((0.93968+0.0622*(x(i)-2.1))/(gam-1))+0.5*gam*((U2(i)/U1(i))^2));
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
        
        dU1dt_p(i)=-(F11-F1)/dx;
        dU2dt_p(i)=-((F21-F2)/dx)+J2;
        dU3dt_p(i)=-(F31-F3)/dx;

        % Artificial viscosity
        S1=0;
        S2=0;
        S3=0;
        if i>1
            p_1=(gam-1)*(U3(i-1)-0.5*gam*((U2(i-1)^2)/U1(i-1)))/A(i-1);
            p1=(gam-1)*(U3(i)-0.5*gam*((U2(i)^2)/U1(i)))/A(i);
            p11=(gam-1)*(U3(i+1)-0.5*gam*((U2(i+1)^2)/U1(i+1)))/A(i+1);

            lambda=Cx*abs(p_1-2*p1+p11)/(p_1+2*p1+p11);

            S1=lambda*(U1(i-1)-2*U1(i)+U1(i+1));
            S2=lambda*(U2(i-1)-2*U2(i)+U2(i+1));
            S3=lambda*(U3(i-1)-2*U3(i)+U3(i+1));
        end

        U1_m(i)=U1(i)+dU1dt_p(i)*dt+S1;
        U2_m(i)=U2(i)+dU2dt_p(i)*dt+S2;
        U3_m(i)=U3(i)+dU3dt_p(i)*dt+S3;
    end

    for i=2:N-1
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

        dU1dt_av=(dU1dt_p(i)+dU1dt_c)/2;
        dU2dt_av=(dU2dt_p(i)+dU2dt_c)/2;
        dU3dt_av=(dU3dt_p(i)+dU3dt_c)/2;

        % Artificial viscosity
        S1=0;
        S2=0;
        S3=0;
        if i<N-1
            p_1=(gam-1)*(U3_m(i-1)-0.5*gam*((U2_m(i-1)^2)/U1_m(i-1)))/A(i-1);
            p1=(gam-1)*(U3_m(i)-0.5*gam*((U2_m(i)^2)/U1_m(i)))/A(i);
            p11=(gam-1)*(U3_m(i+1)-0.5*gam*((U2_m(i+1)^2)/U1_m(i+1)))/A(i+1);

            lambda=Cx*abs(p_1-2*p1+p11)/(p_1+2*p1+p11);

            S1=lambda*(U1_m(i-1)-2*U1_m(i)+U1_m(i+1));
            S2=lambda*(U2_m(i-1)-2*U2_m(i)+U2_m(i+1));
            S3=lambda*(U3_m(i-1)-2*U3_m(i)+U3_m(i+1));
        end

        U1_n1(i)=U1(i)+dU1dt_av*dt+S1;
        U2_n1(i)=U2(i)+dU2dt_av*dt+S2;
        U3_n1(i)=U3(i)+dU3dt_av*dt+S3;

        % Residuals' update
        if(dU1dt_av>r_U1(n))
            r_U1(n)=dU1dt_av;
        end
        if(dU2dt_av>r_U2(n))
            r_U2(n)=dU2dt_av;
        end
        if(dU3dt_av>r_U3(n))
            r_U3(n)=dU3dt_av;
        end
    end
   
    % Subsonic inflow B.Cs.
    U1_n1(1)=A(1);
    U2_n1(1)=2*U2_n1(2)-U2_n1(3);
    U3_n1(1)=U1_n1(1)*((1/(gam-1))+0.5*gam*((U2_n1(1)/U1_n1(1))^2));
    
    % Supersonic outflow B.Cs.
    U1_n1(N)=2*U1_n1(N-1)-U1_n1(N-2);
    U2_n1(N)=2*U2_n1(N-1)-U2_n1(N-2);
    U3_n1(N)=(pepo*A(N)/(gam-1))+0.5*gam*((U2_n1(N)^2)/U1_n1(N));

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

figure(2);
set(gcf,'color','white');
plot(x,Comparator(:,6));

figure(3);
set(gcf,'color','white');
plot(x,Comparator(:,7));

figure(4);
set(gcf,'color','white');
plot(x,Comparator(:,8));
