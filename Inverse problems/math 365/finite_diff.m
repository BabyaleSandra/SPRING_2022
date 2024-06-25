%% Rocket problem
clear;close all;
n=10;
a=0;b=5;
t=linspace(a,b,n);
h=(b-a)/n;
g=-9.8;

u=ones(n,1);
l=u;
u(1:2)=0;
l(end-1:end)=0;
d=-2*ones(n,1);
d(1)=1;
d(end)=1;
A=spdiags([l d u],-1:1,n,n);
%A(1,2)=0;A(end,end)=1;
b=g*h^2*ones(n,1);
b(1)=0;
b(end)=50;
y_est=A\b;

h=plot(t,y_est,'linewidth',2);
h_gca=gca;
h_gca.FontSize=14;
xlabel('time (s)')
ylabel('altitude (m)')
title('Rocket problem')

%% Boundary value problem
clear;close all;
n=6;
a=0;b=pi/2;
x=linspace(a,b,n);
h=(b-a)/n;

u=ones(n,1);
l=u;
u(1:2)=0;
l(end-1)=2;
l(end)=0;
d=(-2+4*h^2)*ones(n,1);
d(1)=1;
A=spdiags([l d u],-1:1,n,n);
% b=zeros(n,1);
% b(2:end)=4*h^2*x(2:end);
b=4*h^2*x';
b(1)=0;
y_est=A\b;

h=plot(x,y_est,'linewidth',2);
h_gca=gca;
h_gca.FontSize=14;
xlabel('x')
ylabel('y')
title('Boundary value problem')
