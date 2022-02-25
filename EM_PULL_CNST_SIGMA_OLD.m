function [ efficiencyScore,allocationScore,reliabilityScore, maxSINR, numberOfServedUsers,iter] = EM_PULL_CNST_SIGMA_OLD (x,m,sigma_in,softmaxScale,maxUsersPerCluster,maxIter,sinr_threshold)
%Here the contribution of x to parameter update of drone j is determined by
ThresholdEm = 0.1;
h=70;
% ThresholdDel = 0.5;
dMuMax = 1000000;
AxisGranularity=2; % number of axis configurations to be tested
n=size(x,1); %number of test points
m=sum(sigma_in~=0);
if m ==0
    allocationScore =100;
    reliabilityScore =0;
    efficiencyScore = 0;
    return
end
mu = x(chooseAxis(x,m,AxisGranularity),:);% this function return the indices of points to be used as initial means
sigma_in(sigma_in==0) = [];
sigma=sigma_in;
% for i =1:length(sigma_in)
%     initialSigma = sigma_in(i);
%     sigma = cat(3,sigma,[initialSigma initialSigma].*eye(2));
% end
% SINR_desired_threshold = db2mag(20); %linear
SINR_EM_Threshold = 0.5;
% maxUsersPerCluster = 500;
deleted = 20;
iter = 0;
% maxIter = 500;
muNew = mu;
softmaxScale = 3;
while iter < maxIter
    iter= iter + 1;
    muT=mu;
    %     input('')
    %     hold off;
    %     plotNetwork(x,sigma./10,mu,[min(x(:,1)) max(x(:,1))],[min(x(:,2)) max(x(:,2))]);
    SINR = getSINRFinal(x,mu,sigma./1000,h.*ones(1,m));
    %     SINR = getSINR(m,x,mu,sigma);
    [maxSINR, indices] = max(SINR,[],2);
    %     usersPerCluster = sum(indices==(1:m));
    T = SINR >= (max(SINR,[],2)- SINR_EM_Threshold);
    T = softmax((SINR)'.*softmaxScale)'.*T;
    T= T./sum(T,2); %normalize posterior probabilities
    Pj =sum(T,1);
    Pj(Pj==0)=1;
    for j=1:m
        muNew(j,:)=[T(:,j)'*x(:,1) T(:,j)'*x(:,2)]./Pj(j);
        dMu = muNew(j,:) - mu(j,:);
        if norm(dMu) > dMuMax
            mu(j,:) = mu(j,:) + dMuMax * dMu./norm(dMu);
        else
            mu(j,:) = muNew(j,:);
        end
    end
    Pj = Pj./n;
    muT= abs(mu-muT);
    maxChange=max(max(muT,[],'all'));
    if maxChange<ThresholdEm
        if deleted > 0
            deleted = deleted -1;
        end
        if deleted==0
            break
        end
    end
end
allocationScore = sum(indices==(1:m));
% numberOfServedUsers = n - sum(allocationScore(sum(indices==(1:m))>maxUsersPerCluster)-maxUsersPerCluster);
usersExceeding = sum(((maxSINR > 10.^(sinr_threshold./10)).*indices)==(1:m))...
    - maxUsersPerCluster;
usersExceeding(usersExceeding<0)=0;
numberOfServedUsers = (sum(maxSINR > 10.^(sinr_threshold./10))- sum(usersExceeding));

allocationScore = mean(allocationScore(sum(indices==(1:m))>maxUsersPerCluster)-maxUsersPerCluster);
if isnan(allocationScore)
    allocationScore = 0;
end
reliabilityScore = sum(maxSINR > 10.^(sinr_threshold./10))/n + 1;
efficiencyScore = (sum(log2(maxSINR + 1)))./sum(sigma_in);