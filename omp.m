%%%%%%%%%%%%%%%%%%%%%%%%% Orthogonal Matching Pursuit
%%%%%%%%%%%%%%%%%%%%%%%%% Inputs: 
% A - measurement matrix
% y - measurements
% sparsity - sparse level of the channel in angle domain
% thres - threshold for OMP
%%%%%%%%%%%%%%%%%%%%%%%%% Outputs: 
% x_est - estimated sparse channel vector
function x_est = omp(A, y, sparsity, thres)


argmtx = [];
indx = [];
res = y;
x_est = zeros(length(A(1,:)),1);

for k = 1:sparsity
    vec = A'*res;
    [~,indx(k)] = max(abs(vec));%maximum projection
    argmtx(:,k) = A(:,indx(k));%basis matrix
    x_val = inv(argmtx' * argmtx) * (argmtx') * y;

    res = y-argmtx * x_val;%residue(update y)
    
    %reconstruct x
    for m = 1: length(indx)
        x_est(indx(m)) = x_val(m);
    end
    
    %criterian
    if norm(res) < thres
        break
    end
end