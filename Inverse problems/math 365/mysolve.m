
A=[1 2 3;2 1 4;3 4 1];
b=[1;1;1];
Ainv=inv(A);
[L,U]=lu(A);
xinv=Ainv*b;
y=L\b; xLU=U\y;
xback=A\b;
% calculate norms
xtrue=[0;0.2;0.2]
xest=[1e-4;0.2;0.19]
norm(xtrue-xest,1)
norm(xtrue-xest,2)
norm(xtrue-xest,inf)

%%
A=[1 0 2;2 -1 3;4 1 8]
Ainv=inv(A);
[L,U,P]=lu(A);
% Identify if LU=A, what is L
%%
n=4000;
S=rand(n,n);
f=rand(n,1);

ntry=50;

tic
for i=1:ntry
    x1=S\f;
end
toc

tic
for i=1:ntry
    x2=inv(S)*f;
end
toc

tic
for i=1:ntry
    [LS,US]=lu(S);
    x3=US\(LS\f);
end
toc
%%
norm(x1-x2)
norm(x1-x3)
norm(x2-x3)
norm(x1-x2,1)
norm(x1-x3,1)
norm(x2-x3,1)
norm(x1-x2,inf)
norm(x1-x3,inf)
norm(x2-x3,inf)



