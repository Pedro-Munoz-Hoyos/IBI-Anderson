% ==================================
% FLOW PARAMETERS
% ==================================
u_fs=1;     % Free-stream velocity [m/s]
R_c=1;      % Radius of the cylinder [m]
R_fs=20;    % Radius of the external domain boundary [m] (R_fs>R_c);

N_r=60;    % Number of radial grid points
N_tht=120;  % Number of angular grid points
w=1.8;      % Relaxation factor;
res=10^(-6);% Convergence criterion over the residuals


% ==================================
% GRID GENERATION
% ==================================
dtht=2*pi/(N_tht-1);
dr=(R_fs-R_c)/(N_r-1);
tht=0:dtht:2*pi;
r=R_c:dr:R_fs;

x=zeros(N_tht,N_r);
y=zeros(N_tht,N_r);
for i=1:N_tht
    for j=1:N_r
        x(i,j)=r(j)*cos(tht(i));
        y(i,j)=r(j)*sin(tht(i));
    end
end

figure(1);
[thtm,rm]=meshgrid(tht,r);
scatter(thtm,rm);
xlabel("\theta [rad]");
ylabel("r [m]");
title("Computational space");
xlim([0 2*pi]);
ylim([R_c R_fs]);

figure(2);
scatter(x,y);
xlabel("x [m]");
ylabel("y [m]");
title("Physical space");
xlim([-R_fs R_fs]);
ylim([-R_fs R_fs]);


% ==================================
% NUMERICAL ENGINE
% ==================================
Phi_n=zeros(N_tht,N_r);
Phi_n1=zeros(N_tht,N_r);

% Upper boundary condition (Dirichlet B.C.)  
for i=1:N_tht
    Phi_n1(i,N_r)=u_fs*R_fs*cos(tht(i));
end

n=1;
nv=[];
mdphiv=[];
maxdphi=res+1;
while maxdphi>res
    maxdphi=0;
    Phi_n=Phi_n1;
    % Lower boundary condition (Neumann B.C.) 
    for i=1:N_tht
        Phi_n(i,1)=(4*Phi_n(i,2)-Phi_n(i,3))/3;
    
        % Check convergence
        if(abs(Phi_n(i,1)-Phi_n1(i,1))>maxdphi)
            maxdphi=abs(Phi_n(i,1)-Phi_n1(i,1));
        end
    end
    
    Phi_n1=Phi_n;
    for i=1:N_tht
        for j=2:N_r-1
            A=(Phi_n1(i,j+1)-Phi_n1(i,j-1))/(2*dr*r(j));
            
            % Lateral boundary condition (0,r)<->(2pi,r) (Neumann B.Cs.) 
            if(i==1)
                B=(Phi_n1(i+1,j)+Phi_n1(N_tht-1,j))/((dtht*r(j))^2);
            elseif(i==N_tht)
                B=(Phi_n1(2,j)+Phi_n1(i-1,j))/((dtht*r(j))^2);
            else
                B=(Phi_n1(i+1,j)+Phi_n1(i-1,j))/((dtht*r(j))^2);
            end
    
            C=(Phi_n1(i,j+1)+Phi_n1(i,j-1))/(dr^2);
            D=2*(((1/(r(j)*dtht))^2)+((1/dr)^2));
            
            Phi_mean=(A+B+C)/D;
            Phi_n1(i,j)=Phi_n(i,j)+w*(Phi_mean-Phi_n(i,j));
    
            % Check convergence
            if(abs(Phi_n1(i,j)-Phi_n(i,j))>maxdphi)
                maxdphi=abs(Phi_n1(i,j)-Phi_n(i,j));
            end
        end
    end
    nv(n)=n;
    mdphiv(n)=maxdphi;
    
    fig3=figure(3);
    refresh(fig3);
    set(gcf,'color','white');
    plot(nv,mdphiv);
    xlabel("Iterations");
    ylabel("max(\Delta\phi)");
    title("Simulation convergence");
    
    fig4=figure(4);
    refresh(fig4);
    set(gcf,'color','white');
    surf(tht,r,transpose(Phi_n));
    xlabel("\theta [rad]");
    ylabel("r [m]");
    zlabel("\phi [m^2/s]");
    title("Velocity potential \phi(r,\theta)");
    xlim([0 2*pi]);
    ylim([R_c R_fs]);
    pause(1/10);
    
    n=n+1;
end


% ==================================
% POST-PROCESSING
% ==================================
u=zeros(N_tht,N_r);
v=zeros(N_tht,N_r);
Cp_n=zeros(N_tht,1);
Cp_a=zeros(N_tht,1);

maxError=0;
for i=1:N_tht
    for j=1:N_r
        Phi_a_ij=u_fs*r(j)*(1+((R_c/r(j))^2))*cos(tht(i)); % Analytical solution
        
        % Obtains the maximum error (absolute) between the analytical and numerical
        % solutions
        if(abs(Phi_n(i,j)-Phi_a_ij)>maxError)
            maxError=(abs(Phi_n(i,j)-Phi_a_ij));
        end

        if(i==1)
            dphi_dtht=(Phi_n(i+1,j)-Phi_n(N_tht-1,j))/(2*dtht);
        elseif(i==N_tht)
            dphi_dtht=(Phi_n(2,j)-Phi_n(i-1,j))/(2*dtht);
        else
            dphi_dtht=(Phi_n(i+1,j)-Phi_n(i-1,j))/(2*dtht);
        end

        if(j==1)
            dphi_dr=(-3*Phi_n(i,j)+4*Phi_n(i,j+1)-Phi_n(i,j+2))/(2*dr);
        elseif(j==N_r)
            dphi_dr=(3*Phi_n(i,j)-4*Phi_n(i,j-1)+Phi_n(i,j-2))/(2*dr);
        else
            dphi_dr=(Phi_n(i,j+1)-Phi_n(i,j-1))/(2*dr);
        end

        u(i,j)=dphi_dtht*(-sin(tht(i))/r(j))+dphi_dr*cos(tht(i));
        v(i,j)=dphi_dtht*(cos(tht(i))/r(j))+dphi_dr*sin(tht(i));
    end

    Cp_n(i)=1-(((u(i,1)^2)+(v(i,1)^2))/(u_fs^2)); % Numerical
    Cp_a(i)=1-((-2*sin(tht(i)))^2);               % Analytical
end

figure(5);
set(gcf,'color','white');
contourf(x,y,u,40);
xlabel("x [m]");
ylabel("y [m]");
title("X-velocity field u[m/s]");
xlim([-R_fs R_fs]);
ylim([-R_fs R_fs]);

figure(6);
set(gcf,'color','white');
plot(tht,Cp_n,'r');
hold on;
plot(tht,Cp_a,'b');
hold off;
xlabel("\theta [rad]");
ylabel("C_p");
title("Pressure coefficient's distribution C_p(\theta)");
xlim([0 2*pi]);
