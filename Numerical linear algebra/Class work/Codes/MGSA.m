% ==============================================================
% This function implement the modified Gram-Schmidt algorithm
%
% Input: 
%       Matrix A
%
% Outputs: 
%       Matrix Q
%       Matrix R, an upper triangular matrix
% ==============================================================
function [Q, R] = MGSA(A)

[m,n] = size(A);

Q = zeros(m,n);
R = zeros(n,n);
v = zeros(m,n);

for j = 1:n
    
    v(:,j) = A(:,j);

end
for i = 1:n

    R(i,i) = norm(v(:,i),2);
    Q(:,i) = v(:,i)/R(i,i);
    
    for j = i+1:n
        
        R(i,j) = Q(:,i)'*v(:,j);
        v(:,j) = v(:,j) - R(i,j)*Q(:,i);
    end
end
end