clear all
close all
clc
% Dimension on the matrix
m = 100;   
n = 15;

% Vendermond matrix
A = zeros(m,n);
I = eye(n);
J = 1:m;
t = (J-1)/(m-1);

for j = 1:n
    
    A(:,j) = t.^(j-1);
end

% call of CGSA
[Q_CG,R_CG] = CGSA(A);
% computation of matrix A = QR using CGSA
A_CG = Q_CG*R_CG;
fprintf('The norms ||A-QR|| and ||I-Q*Q|| in case of classical Gram-Schmidt:\n\n')
% infinty norm for the matrix A and the computed one
Norm_ACG = norm(A-A_CG,inf);
fprintf('||A-QR|| = ')
disp(Norm_ACG)
% infinty norm for the matrix I and Q'*Q
Norm_QCG = norm(I-Q_CG'*Q_CG, inf);
fprintf('||I-Q*Q|| = ')
disp(Norm_QCG)
% call of MGSA
[Q_MG,R_MG] = MGSA(A);
% computation of matrix A = QR using MGSA
A_MG = Q_MG*R_MG;
fprintf('The norms ||A-QR|| and ||I-Q*Q|| in case of modified Gram-Schmidt:\n\n')
% infinty norm for the matrix A and the computed one
Norm_AMG = norm(A-A_MG,inf);
fprintf('||A-QR|| =')
disp(Norm_AMG)
% infinty norm for the matrix I and Q'*Q
Norm_QMG = norm(I-Q_MG'*Q_MG, inf);
fprintf('||I-Q*Q|| =')
disp(Norm_QMG)

disp('From above results one can see that ||A-QR|| is almoast zero for both')
disp('algorithms. But the strang thing is that the norm ||I-Q*Q|| for ')
disp('classical algorithm is much bigger while for the modified algorithm it')
disp('is very small. So one can conclude that the classical method is unstable.')