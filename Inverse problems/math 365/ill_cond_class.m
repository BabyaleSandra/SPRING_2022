% Examples
clear;close all

% Grams of three different compounds A,B,C held in a single solution.
disp('-------------------------')
disp('Compound example')
disp('-------------------------')
pause
A=[93.477 10.202 -28.832;1.963 32.816 62.414;26.821 36.816 57.234]
pause
disp(['Matrix rank is ', num2str(rank(A))]);
pause
b=[34.7177;70.9241;82.9271]
pause
x=A\b
pause
b_new=[34.7;70.9;82.9]
pause
x_new=A\b_new
pause
disp('                                  ')
%%
% Golub examples
disp('-------------------------')
disp('Golub example')
disp('-------------------------')
pause
prompt='Choose a value for n  '
n=input(prompt);
AG=golub(n);
disp(['Matrix rank is ', num2str(rank(AG))]);
pause
xG=rand(n,1);
bG=AG*xG;
xG_est=AG\bG;
pause
error=norm(xG_est-xG);
disp(['|| error ||_2 = ',num2str(error)]);
pause
disp('                                  ')
%% Error in compound example
disp('-------------------------')
disp('Compound example')
disp('-------------------------')
pause
delb=b_new-b;
rel_b=norm(delb)/norm(b);
disp(['Relative error in b = ',num2str(rel_b)]);
pause
delx=x_new-x;
rel_x=norm(delx)/norm(x);
disp(['Relative error in x = ',num2str(rel_x)]);
pause
min_cond=rel_x/rel_b;
disp(['Lower bound on condition number of A = ',num2str(min_cond)]);
pause
cond(A);
disp(['Condition number of A = ',num2str(cond(A))]);
pause
disp('                                  ')
%% Error in Golub example
disp('-------------------------')
disp('Golub example')
disp('-------------------------')
pause
bG_new=bG+0.02*rand(size(bG));
delbG=bG_new-bG;
rel_bG=norm(delbG)/norm(bG);
disp(['Relative error in b = ',num2str(rel_bG)]);
pause
xG_new=AG\bG_new;
pause
delxG=xG_new-xG;
rel_xG=norm(delxG)/norm(xG);
disp(['Relative error in x = ',num2str(rel_xG)]);
pause
min_cond=rel_xG/rel_bG;
disp(['Lower bound on condition number of A = ',num2str(min_cond)]);
pause
cond(AG);
disp(['Condition number of A = ',num2str(cond(AG))]);
