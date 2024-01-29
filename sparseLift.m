%%%%%%%%%%%%%%%%%%%%%%%%% SparseLift
%%%%%%%%%%%%%%%%%%%%%%%%% Inputs: 
% A - measurement matrix
% B - structured matrix(diagnal matrix of all-1-vectors)
% y - measurements
% sparsity - sparse level of the channel in angle domain
% thres - threshold for OMP
%%%%%%%%%%%%%%%%%%%%%%%%% Outputs: 
% x_est - estimated sparse channel vector
function x_est = sparseLift(A,B,y,pack_num,sparse,thres)
m = length(A(:,1));
n = length(A(1,:));
x_est = zeros(n,1);

phi_operation = zeros(pack_num*n,m);
for k = 1:1:m
    phi_operation(:,k) = kron((A(k,:))',(B(k,:))'); 
end
clear k;
phi_operation = phi_operation';
% reconstruct
z = omp(phi_operation,y,sparse*pack_num,thres);
% v_est = ISTA(phi_operation,y,z,sparse*pack_num,1);

% estimate x
X_est = reshape(z, pack_num, n);
[~,S,V] = svd(X_est);
x_est_pre = sqrt(S(1,1))* (V(:,1))';
indices=retin(abs(x_est_pre),sparse);
x_est(indices) = x_est_pre(indices);
end
