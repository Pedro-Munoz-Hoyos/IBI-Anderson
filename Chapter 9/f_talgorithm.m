function [un1,res] = talgorithm(u,E,dt)
% Function in charge of implementing Thomas' algorithm.

N=length(u);
A=-0.5*E;
B=1+E;

un1=zeros(N,1);
Dp=zeros(N-2,1);
Cp=zeros(N-2,1);

% Transformation of the tridiagonal form into upper bidiagonal.
Dp(1)=B;
Cp(1)=((1-E)*u(2))-A*(u(1)+u(3))-A*u(1);
for j=2:N-2
    Dp(j)=B-((A^2)/Dp(j-1));
    Cp(j)=(1-E)*u(j+1)-A*(u(j)+u(j+2))-(Cp(j-1)*A/Dp(j-1));
end
Cp(N-2)=Cp(N-2)-A*u(N);

% Upward resolution of the upper bidiagonal system.
un1(N)=u(N);
un1(N-1)=Cp(N-2)/Dp(N-2);
res=abs(un1(N-1)-u(N-1));
for j=2:N-2
    un1(N-j)=(Cp(N-j-1)-A*un1(N-j+1))/Dp(N-j-1);
    dudt=abs(un1(N-j)-u(N-j));
    if(dudt>res)
        res=dudt; % Residuals' update
    end
end
res=res/dt;

end
