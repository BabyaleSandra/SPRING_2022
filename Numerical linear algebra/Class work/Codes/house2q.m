%======================================================================
% This function takes as input the matrix V from Householder QR
% decomposition and m-by-m orthogonal matrix Q
%======================================================================

function Q = house2q(V)

[m,n] = size(V);

% Initialization of Q to an identity matrix
Q = eye(m);

% Computation of Q using implicit calculation of the product Qx
% x is represented in this context by the columns of the identity matrix

for k = n:-1:1

    Q(k:m,:) = Q(k:m,:) - 2*V(k:m,k)*(V(k:m,k)'*Q(k:m,:));

end
end