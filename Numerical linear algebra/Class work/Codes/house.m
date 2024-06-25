%======================================================================
% This function is built base on Householder QR decomposition algorithm
% Its takes as inputs, matrix A and return matrix V and R
%======================================================================


function [V,R] = house(A)

[m,n] = size(A);

V = zeros(m,n);

for k = 1:n
    
    x = A(k:m,k);
    % The computation of vk = sign(x(1))*norm(x,2)*eye(m-k+1,1) + x;
    % can be done by adding "sign(x(1))*norm(x,2)" to x(1) the first
    % element of x and consider vk as the new x.
    
    x(1) = x(1)+  sign(x(1))*norm(x,2);
    vk = x;
    vk = vk/norm(vk,2);
    V(k:m,k) = vk;
    
    A(k:m,k:n) = A(k:m,k:n) - 2*vk*(vk'*A(k:m,k:n));
        
end

R = A;
end