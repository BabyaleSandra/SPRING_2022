n=51;
W=145; S=27;xin=0.15; yin=0; m=7;
u=S*m*ones(n,1);
l=W*ones(n,1);
u(1)=0;
l(end)=0;
d=-(W+S*m)*ones(n,1);
A=spdiags([l d u],-1:1,n,n);
b=zeros(n,1);
b(1)=-W*xin;
b(n)=-S*yin;
x=A\b;
y=m*x;
plot(x);hold on;plot(y);legend('x','y')
% Amount change
abs(xin-x(end))/xin
