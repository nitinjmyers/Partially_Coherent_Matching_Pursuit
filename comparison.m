clear all;

simulation_id=3; % 1- Performance with SNR, 2- Performance with CS measurements, 3- Performance with sparsity level 

if(simulation_id==1)
    %% white noise power changes from -10dB~15dB
    packetSize = 16;
    sparsity = 4; %sparse level
    n = 256; %N antennas
    m = 128; %M measurements
    packetNumber = m/packetSize;
    aveNumber = 40; % Execution times for averaging (Used 400 for the ICC paper)
    fc60GHz = 60e+9; %60Ghz
    jiterConst = 4.7e-18; %jiter constant
    
    count = 1;
    for snr = -10:5:15
        addNoiPow = 10^(-snr/10); %power of additive noise
        thres = sqrt(addNoiPow);
        error_std = zeros(aveNumber,1); error_omp = zeros(aveNumber,1); error_lift = zeros(aveNumber,1);
        error_pccpr = zeros(aveNumber,1); error_pcmp = zeros(aveNumber,1); error_pccpr10 = zeros(aveNumber,1);
        achievRate_std = zeros(aveNumber,1);    achievRate_omp = zeros(aveNumber,1);    achievRate_lift = zeros(aveNumber,1);
        achievRate_pcmp = zeros(aveNumber,1);    achievRate_pccpr = zeros(aveNumber,1);    achievRate_pccpr10 = zeros(aveNumber,1);
        for t = 1:1:aveNumber
            [AA,A,B,yy,yy_mis,y,y_mis,y_nonPhase,A_p,y_p,y_p_mis,x] = data_gen(m,n,sparsity,packetSize,fc60GHz,jiterConst,addNoiPow);
            x_std = omp(AA,y_nonPhase,sparsity,thres);
            x_omp = omp(AA,yy_mis,sparsity,thres);
            x_lift = sparseLift(AA,B,yy_mis,packetNumber,sparsity,thres);
            x_pcmp = PC_MP(A,A_p,y_p_mis,sparsity,packetNumber);
            x_pccpr = PC_CPR(A,y_mis,A_p,y_p_mis,sparsity,packetNumber,thres,1);
            x_pccpr10 = PC_CPR(A,y_mis,A_p,y_p_mis,sparsity,packetNumber,thres,10);
            %normalized
            x_std = x_std/norm(x_std);
            x_omp = x_omp/norm(x_omp);
            x_lift = x_lift/norm(x_lift);
            x_pcmp = x_pcmp/norm(x_pcmp);
            x_pccpr = x_pccpr/norm(x_pccpr);
            x_pccpr10 = x_pccpr10/norm(x_pccpr10);
            
            pha_omp_conpen = exp(1j*angle(x'*x_omp));
            pha_lift_conpen = exp(1j*angle(x'*x_lift));
            pha_pcmp_conpen = exp(1j*angle(x'*x_pcmp));
            pha_pccpr_conpen = exp(1j*angle(x'*x_pccpr));
            pha_pccpr_conpen10 = exp(1j*angle(x'*x_pccpr10));
            error_std(t) = 10*log10(norm(x-x_std)/norm(x));
            error_omp(t) = 10*log10(norm(x*(pha_omp_conpen)-x_omp)/norm(x));
            error_lift(t) = 10*log10(norm(x*(pha_lift_conpen)-x_lift)/norm(x));
            error_pcmp(t) = 10*log10(norm(x*(pha_pcmp_conpen)-x_pcmp)/norm(x));
            error_pccpr(t) = 10*log10(norm(x*(pha_pccpr_conpen)-x_pccpr)/norm(x));
            error_pccpr10(t) = 10*log10(norm(x*(pha_pccpr_conpen10)-x_pccpr10)/norm(x));
            achievRate_std(t) = achiRateCalc(x,x_std,addNoiPow);
            achievRate_omp(t) = achiRateCalc(x,x_omp,addNoiPow);
            achievRate_lift(t) = achiRateCalc(x,x_lift,addNoiPow);
            achievRate_pcmp(t) = achiRateCalc(x,x_pcmp,addNoiPow);
            achievRate_pccpr(t) = achiRateCalc(x,x_pccpr,addNoiPow);
            achievRate_pccpr10(t) = achiRateCalc(x,x_pccpr10,addNoiPow);
        end
        clear t
        err_pow_std(count) = sum(error_std)/aveNumber;
        err_pow_omp(count) = sum(error_omp)/aveNumber;
        err_pow_lift(count) = sum(error_lift)/aveNumber;
        err_pow_pcmp(count) = sum(error_pcmp)/aveNumber;
        err_pow_pccpr(count) = sum(error_pccpr)/aveNumber;
        err_pow_pccpr10(count) = sum(error_pccpr10)/aveNumber;
        R_pow_std(count) = sum(achievRate_std)/aveNumber;
        R_pow_omp(count) = sum(achievRate_omp)/aveNumber;
        R_pow_lift(count) = sum(achievRate_lift)/aveNumber;
        R_pow_pcmp(count) = sum(achievRate_pcmp)/aveNumber;
        R_pow_pccpr(count) = sum(achievRate_pccpr)/aveNumber;
        R_pow_pccpr10(count) = sum(achievRate_pccpr10)/aveNumber;
        
        count = count + 1;
    end
    
end

if(simulation_id==2)
    %% measurements changes from 16~144
    snr = 15;
    addNoiPow = 10^(-snr/10);
    packetSize = 16;
    sparsity = 4;
    n = 256;
    fc60GHz = 60e+9;
    jiterConst = 4.7e-18;
    aveNumber = 40; % Used 400 for the ICC paper
    thres = 0.1;
    
    err_mes_std = []; err_mes_omp = []; err_mes_lift = []; err_mes_pcmp = []; err_mes_pccpr = []; err_mes_pccpr10 = [];
    count = 1;
    for m = 16:32:144
        packetNumber = m/packetSize;
        error_std = zeros(aveNumber,1); error_omp = zeros(aveNumber,1); error_lift = zeros(aveNumber,1);
        error_pccpr = zeros(aveNumber,1); error_pcmp = zeros(aveNumber,1); error_pccpr10 = zeros(aveNumber,1);
        achievRate_std = zeros(aveNumber,1);    achievRate_omp = zeros(aveNumber,1);    achievRate_lift = zeros(aveNumber,1);
        achievRate_pcmp = zeros(aveNumber,1);    achievRate_pccpr = zeros(aveNumber,1);    achievRate_pccpr10 = zeros(aveNumber,1);
        for t = 1:1:aveNumber
            [AA,A,B,yy,yy_mis,y,y_mis,y_nonPhase,A_p,y_p,y_p_mis,x] = data_gen(m,n,sparsity,packetSize,fc60GHz,jiterConst,addNoiPow);
            x_std = omp(AA,y_nonPhase,sparsity,thres);
            x_omp = omp(AA,yy_mis,sparsity, thres);
            x_lift = sparseLift(AA,B,yy_mis,packetNumber,sparsity,thres);
            x_pcmp = PC_MP(A,A_p,y_p_mis,sparsity,packetNumber);
            x_pccpr = PC_CPR(A,y_mis,A_p,y_p_mis,sparsity,packetNumber,thres,1);
            x_pccpr10 = PC_CPR(A,y_mis,A_p,y_p_mis,sparsity,packetNumber,thres,10);
            %normalized
            x_std = x_std/norm(x_std);
            x_omp = x_omp/norm(x_omp);
            x_lift = x_lift/norm(x_lift);
            x_pcmp = x_pcmp/norm(x_pcmp);
            x_pccpr = x_pccpr/norm(x_pccpr);
            x_pccpr10 = x_pccpr10/norm(x_pccpr10);
            pha_omp_conpen = (x_omp'*x_omp)\x_omp' * x;
            pha_lift_conpen = (x_lift'*x_lift)\x_lift' * x;
            pha_pcmp_conpen = (x_pcmp'*x_pcmp)\x_pcmp' * x;
            pha_pccpr_conpen = (x_pccpr'*x_pccpr)\x_pccpr' * x;
            pha_pccpr10_conpen = (x_pccpr10'*x_pccpr10)\x_pccpr10' * x;
            error_std(t) = 10*log10(norm(x-x_std)/norm(x));
            error_omp(t) = 10*log10(norm(x-x_omp*(pha_omp_conpen))/norm(x));
            error_lift(t) = 10*log10(norm(x-x_lift*(pha_lift_conpen))/norm(x));
            error_pcmp(t) = 10*log10(norm(x-x_pcmp*(pha_pcmp_conpen))/norm(x));
            error_pccpr(t) = 10*log10(norm(x-x_pccpr*(pha_pccpr_conpen))/norm(x));
            error_pccpr10(t) = 10*log10(norm(x-x_pccpr10*(pha_pccpr10_conpen))/norm(x));
            achievRate_std(t) = achiRateCalc(x,x_std,addNoiPow);
            achievRate_omp(t) = achiRateCalc(x,x_omp,addNoiPow);
            achievRate_lift(t) = achiRateCalc(x,x_lift,addNoiPow);
            achievRate_pcmp(t) = achiRateCalc(x,x_pcmp,addNoiPow);
            achievRate_pccpr(t) = achiRateCalc(x,x_pccpr,addNoiPow);
            achievRate_pccpr10(t) = achiRateCalc(x,x_pccpr10,addNoiPow);
        end
        clear t
        err_mes_std(count) = sum(error_std)/aveNumber;
        err_mes_omp(count) = sum(error_omp)/aveNumber;
        err_mes_lift(count) = sum(error_lift)/aveNumber;
        err_mes_pcmp(count) = sum(error_pcmp)/aveNumber;
        err_mes_pccpr(count) = sum(error_pccpr)/aveNumber;
        err_mes_pccpr10(count) = sum(error_pccpr10)/aveNumber;
        R_mes_std(count) = sum(achievRate_std)/aveNumber;
        R_mes_omp(count) = sum(achievRate_omp)/aveNumber;
        R_mes_lift(count) = sum(achievRate_lift)/aveNumber;
        R_mes_pcmp(count) = sum(achievRate_pcmp)/aveNumber;
        R_mes_pccpr(count) = sum(achievRate_pccpr)/aveNumber;
        R_mes_pccpr10(count) = sum(achievRate_pccpr10)/aveNumber;
        
        count = count + 1;
    end
    
end

if(simulation_id==3)
    %% sparse changes from 1~128
    snr = 15;
    addNoiPow = 10^(-snr/10);
    packetSize = 16;
    n = 256;
    m = 128;
    fc60GHz = 60e+9;
    jiterConst = 4.7e-18;
    aveNumber = 40; % Used 400 in the ICC paper
    thres = 0.1;
    
    err_spar_std = []; err_spar_omp = []; err_spar_lift = []; err_spar_pcmp = []; err_spar_pccpr = []; err_spar_pccpr10 = [];
    count = 1;
    for sparsity = [1,2,4,6,8]
        packetNumber = m/packetSize;
        error_std = zeros(aveNumber,1); error_omp = zeros(aveNumber,1); error_lift = zeros(aveNumber,1);
        error_pccpr = zeros(aveNumber,1); error_pcmp = zeros(aveNumber,1); error_pccpr10 = zeros(aveNumber,1);
        achievRate_std = zeros(aveNumber,1);    achievRate_omp = zeros(aveNumber,1);    achievRate_lift = zeros(aveNumber,1);
        achievRate_pcmp = zeros(aveNumber,1);    achievRate_pccpr = zeros(aveNumber,1);    achievRate_pccpr10 = zeros(aveNumber,1);
        for t = 1:1:aveNumber
            [AA,A,B,yy,yy_mis,y,y_mis,y_nonPhase,A_p,y_p,y_p_mis,x] = data_gen(m,n,sparsity,packetSize,fc60GHz,jiterConst,addNoiPow);
            x_std = omp(AA,y_nonPhase,sparsity,thres);
            x_omp = omp(AA,yy_mis,sparsity, thres);
            x_lift = sparseLift(AA,B,yy_mis,packetNumber,sparsity,thres);
            x_pcmp = PC_MP(A,A_p,y_p_mis,sparsity,packetNumber);
            x_pccpr = PC_CPR(A,y_mis,A_p,y_p_mis,sparsity,packetNumber,thres,1);
            x_pccpr10 = PC_CPR(A,y_mis,A_p,y_p_mis,sparsity,packetNumber,thres,10);
            %normalized
            x_std = x_std/norm(x_std);
            x_omp = x_omp/norm(x_omp);
            x_lift = x_lift/norm(x_lift);
            x_pcmp = x_pcmp/norm(x_pcmp);
            x_pccpr = x_pccpr/norm(x_pccpr);
            x_pccpr10 = x_pccpr10/norm(x_pccpr10);
            pha_omp_conpen = (x_omp'*x_omp)\x_omp' * x;
            pha_lift_conpen = (x_lift'*x_lift)\x_lift' * x;
            pha_pcmp_conpen = (x_pcmp'*x_pcmp)\x_pcmp' * x;
            pha_pccpr_conpen = (x_pccpr'*x_pccpr)\x_pccpr' * x;
            pha_pccpr10_conpen = (x_pccpr10'*x_pccpr10)\x_pccpr10' * x;
            error_std(t) = 10*log10(norm(x-x_std)/norm(x));
            error_omp(t) = 10*log10(norm(x-x_omp*(pha_omp_conpen))/norm(x));
            error_lift(t) = 10*log10(norm(x-x_lift*(pha_lift_conpen))/norm(x));
            error_pcmp(t) = 10*log10(norm(x-x_pcmp*(pha_pcmp_conpen))/norm(x));
            error_pccpr(t) = 10*log10(norm(x-x_pccpr*(pha_pccpr_conpen))/norm(x));
            error_pccpr10(t) = 10*log10(norm(x-x_pccpr10*(pha_pccpr10_conpen))/norm(x));
            achievRate_std(t) = achiRateCalc(x,x_std,addNoiPow);
            achievRate_omp(t) = achiRateCalc(x,x_omp,addNoiPow);
            achievRate_lift(t) = achiRateCalc(x,x_lift,addNoiPow);
            achievRate_pcmp(t) = achiRateCalc(x,x_pcmp,addNoiPow);
            achievRate_pccpr(t) = achiRateCalc(x,x_pccpr,addNoiPow);
            achievRate_pccpr10(t) = achiRateCalc(x,x_pccpr10,addNoiPow);
        end
        clear t
        err_spar_std(count) = sum(error_std)/aveNumber;
        err_spar_omp(count) = sum(error_omp)/aveNumber;
        err_spar_lift(count) = sum(error_lift)/aveNumber;
        err_spar_pcmp(count) = sum(error_pcmp)/aveNumber;
        err_spar_pccpr(count) = sum(error_pccpr)/aveNumber;
        err_spar_pccpr10(count) = sum(error_pccpr10)/aveNumber;
        R_spar_std(count) = sum(achievRate_std)/aveNumber;
        R_spar_omp(count) = sum(achievRate_omp)/aveNumber;
        R_spar_lift(count) = sum(achievRate_lift)/aveNumber;
        R_spar_pcmp(count) = sum(achievRate_pcmp)/aveNumber;
        R_spar_pccpr(count) = sum(achievRate_pccpr)/aveNumber;
        R_spar_pccpr10(count) = sum(achievRate_pccpr10)/aveNumber;
        count = count + 1;
    end
end
%%
set(0,'DefaultLineLineWidth',1.5)
set(0,'defaulttextinterpreter','tex')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultLineMarkerSize',6)
set(0,'DefaultAxesFontWeight','normal')
set(gca,'FontSize',12)
set(get(gca,'Xlabel'),'FontSize',10)
set(get(gca,'Ylabel'),'FontSize',10)
%%
if(simulation_id==1)
    figure(1)
    x_axis = -10:5:15;
    plot(x_axis,R_pow_std,'-*')
    hold on
    plot(x_axis,R_pow_pcmp,'-s')
    plot(x_axis,R_pow_omp,'-o')
    plot(x_axis,R_pow_lift,'-+')
    plot(x_axis,R_pow_pccpr,'-d')
    plot(x_axis,R_pow_pccpr10,'-x')
    ylim([0 5.5])
    %xlim([0 5.5])
    legend('OMP-genie phase','PC-MP','OMP','SparseLift','PC-CPR','x10 PC-CPR')
    hold off
    grid on
    xlabel('SNR (dB)')
    ylabel('Achievable Rate (bps/Hz)')
    set(gcf,'unit', 'centimeters','position',[2 5 12 9])
    
    figure(2)
    x_axis = -10:5:15;
    plot(x_axis,err_pow_std,'-*')
    hold on
    plot(x_axis,err_pow_pcmp,'-s')
    plot(x_axis,err_pow_omp,'-o')
    plot(x_axis,err_pow_lift,'-+')
    plot(x_axis,err_pow_pccpr,'-d')
    plot(x_axis,err_pow_pccpr10,'-x')
    grid on
    legend('OMP-genie phase','PC-MP','OMP','SparseLift','PC-CPR','x10 PC-CPR')
    ylim([-23 1.5])
    hold off
    xlabel('SNR (dB)')
    ylabel('NMSE (dB)')
    set(gcf,'unit', 'centimeters','position',[14 5 12 9])
end
%%
if(simulation_id==2)
    figure(1)
    x_axis = 16:32:144;
    plot(x_axis,R_mes_std,'-*')
    hold on
    plot(x_axis,R_mes_pcmp,'-s')
    plot(x_axis,R_mes_omp,'-o')
    plot(x_axis,R_mes_lift,'-+')
    plot(x_axis,R_mes_pccpr,'-d')
    plot(x_axis,R_mes_pccpr10,'-x')
    grid on
    legend('OMP-genie phase','PC-MP','OMP','SparseLift','PC-CPR','x10 PC-CPR')
    xlim([16 144])
    ylim([0.8 5.1])
    hold off
    xlabel('Number of measurements (M·P)')
    ylabel('Achievable Rate (bps/Hz)')
    set(gcf,'unit', 'centimeters','position',[8 5 12 9])
    
    figure(2)
    x_axis = 16:32:144;
    plot(x_axis,err_mes_std,'-*')
    hold on
    plot(x_axis,err_mes_pcmp,'-s')
    plot(x_axis,err_mes_omp,'-o')
    plot(x_axis,err_mes_lift,'-+')
    plot(x_axis,err_mes_pccpr,'-d')
    plot(x_axis,err_mes_pccpr10,'-x')
    grid on
    xlim([16 144])
    ylim([-23 0])
    legend('OMP-genie phase','PC-MP','OMP','SparseLift','PC-CPR','x10 PC-CPR')
    hold off
    xlabel('Number of measurements (M·P)')
    ylabel('NMSE (dB)')
    set(gcf,'unit', 'centimeters','position',[14 5 12 9])
end

if(simulation_id==3)
    figure(1)
    x_axis = [2,4,6,8];
    plot(x_axis,R_spar_std(2:5),'-*')
    hold on
    plot(x_axis,R_spar_pcmp(2:5),'-s')
    plot(x_axis,R_spar_omp(2:5),'-o')
    plot(x_axis,R_spar_lift(2:5),'-+')
    plot(x_axis,R_spar_pccpr(2:5),'-d')
    plot(x_axis,R_spar_pccpr10(2:5),'-x')
    grid on
    legend('OMP-genie phase','PC-MP','OMP','SparseLift','PC-CPR','x10 PC-CPR')
    hold off
    xlabel('Sparsity')
    ylabel('Achievable Rate (bps/Hz)')
    ylim([0.4 5.5])
    set(gcf,'unit', 'centimeters','position',[14 5 12 9])
    
    figure(2)
    x_axis = [2,4,6,8];
    plot(x_axis,err_spar_std(2:5),'-*')
    hold on
    plot(x_axis,err_spar_pcmp(2:5),'-s')
    plot(x_axis,err_spar_omp(2:5),'-o')
    plot(x_axis,err_spar_lift(2:5),'-+')
    plot(x_axis,err_spar_pccpr(2:5),'-d')
    plot(x_axis,err_spar_pccpr10(2:5),'-x')
    grid on
    legend('OMP-genie phase','PC-MP','OMP','SparseLift','PC-CPR','x10 PC-CPR')
    hold off
    xlabel('Sparsity')
    ylabel('NMSE (dB)')
    ylim([-12 0])
    set(gcf,'unit', 'centimeters','position',[14 5 12 9])
end
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%% Partially Coherent Compressive Phase Retrieval
%%%%%%%%%%%%%%%%%%%%%%%%% Inputs:
% A - measurement matrix
% y - measurements
% A_p - measurement matrix corresponding to each packet
% y_p - measurements in each packet
% sparsity - sparse level of the channel in angle domain
% packetNumber - number of the packets
% thres - threshold for hardthresholding
% xtimes - iteration times
%%%%%%%%%%%%%%%%%%%%%%%%% Outputs:
% h_est - estimated sparse channel vector
function h_est = PC_CPR(A,y,A_p,y_p,sparsity,packetNumber,thres,xtimes)
M = size(A,1);
N = size(A,2);
pack_size = M/packetNumber;
% index detection
vec = zeros(N,1);
for b = 1:1:packetNumber
    vec = vec + (abs(A_p(:,:,b)' *y_p(:,b))).^2;
end
[~,idx] = maxk(vec,sparsity);
clear b vec

%initialization
A = A';
pre = zeros(M,1);
r = sqrt(1/M*sum((abs(y)).^2));
for m = 1:1:M
    pre(m) = CRAF(y(m),r);
end
clear m r
Ds = A(idx,:)*diag(pre)*(A(idx,:))';
[d,~]=eigs(Ds,1); %dominant eigenvector
h_ini = (abs((A(idx,:))'*d)).' *abs(y)/((norm((A(idx,:))'*d))^2);
h = zeros(N,1);
h(idx) = h_ini*d;

%refinement
beta = 0.1;
tau = 0.1;
step = 0.2;
h_est = h;
A = A';
cri = norm(y-A*h_est);
t = 0;
while(cri>thres)
    for b = 1:1:packetNumber
        %update y
        y_p(:,b) = y_p(:,b)*(y_p(:,b))'*A_p(:,:,b)*h_est/abs((y_p(:,b))'*A_p(:,:,b)*h_est);
        y(pack_size*(b-1)+1:pack_size*b) = y_p(:,b);
    end
    %update h
    vec_w = abs(A*h)./(abs(A*h)+beta*abs(y));
    vec_w = max(vec_w,tau);
    w = diag(vec_w);
    grad = h_est-(step/M)*(A')*w*(A*h_est-y);
    h_est = hard_thres(grad,sparsity);
    t = t+1;
    if(t>sparsity*xtimes)
        break
    end
end
end

function c = CRAF(ym,r)
lam = 1;% Set lambda = 1 according to the paper
com = (abs(ym))^2;
if com <= r^2/2
    c = -lam/com;
else
    c = lam/com;
end
end

function v = hard_thres(h,sparse)
v = zeros(size(h));
indices=retin(abs(h),sparse);
v(indices) = h(indices);
end

%%%%%%%%%%%%%%%%%%%%%%%%% Achievable Rate with Conjugate Beamformer
function r = achiRateCalc(x,x_est,pow)
beam = x_est';
r = log2(1+abs(beam*x)^2/pow);
end
