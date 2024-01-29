%%%%%%%%%%%%%%%%%%%%%%%%% Partially Coherent Matching Pursuit
%%%%%%%%%%%%%%%%%%%%%%%%% Inputs: 
% A - measurement matrix
% A_p - measurement matrix corresponding to each packet
% y_p - measurements in each packet
% sparsity - sparse level of the channel in angle domain
% packetNumber - number of the packets
%%%%%%%%%%%%%%%%%%%%%%%%% Outputs: 
% x_est - estimated sparse channel vector
function x_est = PC_MP(A,A_p,y_p,sparsity,packetNumber)

n = length(A(1,:));
vec = zeros(n,1);
indx = [];
res = y_p;
x_est = zeros(n,1);
phase_est = randn(packetNumber,1);
for k = 1:sparsity
    %determine index
    for p = 1:1:packetNumber
        vec = vec + abs(A_p(:,:,p)' * res(:,p));
    end
    %different index every time
    [~,idx] = max(vec);
    while (ismember(idx,indx))
        vec(idx) = 0;
        [~,idx] = max(vec);
    end
    indx(k) = idx;
    clear p
    
    %initialize phase simply 
    x_val = zeros(k,1);
    AAmtx = zeros(k,k);
    for p = 1:1:packetNumber
        AAmtx = AAmtx + A_p(:,indx,p)'*A_p(:,indx,p);%can be calculated outside
    end
    %loop for alternating optimization
    for l = 2:1:k*10
        eAy = zeros(k,1);
        %fix x
        for p = 1:1:packetNumber
            phase_est(p) = angle((A_p(:,indx,p)*x_val)' * y_p(:,p));
            phase_est(p) = exp(1j*phase_est(p));
            %fix phase
            eAy = eAy + A_p(:,indx,p)' * y_p(:,p) * (phase_est(p))';
        end
        x_val = (AAmtx'*AAmtx)\(AAmtx') * eAy;%persudo inverse
    end
    %reconstruct x
    for m = 1: length(indx)
        x_est(indx(m)) = x_val(m);
    end
    %update residue for each packet
    for p = 1:1:packetNumber
        res(:,p) = y_p(:,p)-A_p(:,indx,p) * x_val * (phase_est(p)');
    end
    clear p
end
clear k

