% ==============================================================
% This function implement the Classical Gram-Schmidt algorithm
%
% Input: 
%       Matrix A
%
% Outputs: 
%       Matrix Q
%       Matrix R, an upper triangular matrix
% ==============================================================

function [Q, R] = CGSA(A)

[m,n] = size(A);

Q = zeros(m,n);
R = zeros(n,n);

for j = 1:n
    
    vj = A(:,j);
    
    for i = 1:j-1
        
        R(i,j) = Q(:,i)'*A(:,j);
        vj = vj - R(i,j)*Q(:,i);
    
    end
    
    R(j,j) = norm(vj, 2);
    Q(:,j) = vj/R(j,j);   
end
end
   