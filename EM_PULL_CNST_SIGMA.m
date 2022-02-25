function [allocationScore,reliabilityScore] = EM_PULL_CNST_SIGMA (x,m,sigma_in,softmaxScale,maxUsersPerCluster,maxIter)
%Here the contribution of x to parameter update of drone j is determined by
ThresholdEm = 0.1;
ThresholdDel = 0.5;
dMuMax = 1;
AxisGranularity=4; % number of axis configurations to be tested
n=size(x,1); %number of test points
m=sum(sigma_in~=0);
sigma_in(sigma_in==0) = [];
mu = x(chooseAxis(x,m,AxisGranularity),:);% this function return the indices of points to be used as initial means
sigma=[];
for i =1:length(sigma_in)
    initialSigma = sigma_in(i);
    sigma = cat(3,sigma,[initialSigma initialSigma].*eye(2));
end
% initialSigma = sigma;
% sigma = [initialSigma initialSigma].*eye(2); %initialize sigma
% sigma = repmat(sigma,1,1,m);
% sigmaValPrev=initialSigma*ones(1,m);
SINR_desired_threshold = db2mag(20); %linear
% maxUsersPerCluster = 500;
deleted = 20;
iter = 0;
% maxIter = 500;
muNew = mu;
while iter < maxIter
    iter= iter + 1;
    muT=mu;
%     input('')
%     hold off;
%     plotNetwork(x,sigma./100,mu,[min(x(:,1)) max(x(:,1))],[min(x(:,2)) max(x(:,2))]);
    SINR = getSINR(m,x,mu,sigma);
    [maxSINR, indices] = max(SINR,[],2);
    usersPerCluster = sum(indices==(1:m));
    T = logical( maxSINR < SINR_desired_threshold).*SINR + logical(maxSINR > SINR_desired_threshold).*de2bi(indices,m).*maxSINR;
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
    if maxChange<ThresholdDel || iter == maxIter - 40
        ii=0;
        for i=1:m
            iii=i-ii;
%             overlaps = circleOverlapArea(vecnorm(mu(iii,:)-mu,2,2),sqrt(sigma(1,1,iii)./1000),sqrt(squeeze(sigma(1,1,:))./1000));
%             overlaps(isnan(overlaps))=[];
            if any(vecnorm(mu(iii,:)-mu,2,2)<sqrt(squeeze(-sigma(1,1,iii) + sigma(1,1,:))/10)) || sigma(1,1,iii) < 50
                m=m-1;
                mu(iii,:)=[];
                muT(iii,:)=[];
                sigma(:,:,iii)=[];
                Pj(iii)=[];
%                 sigmaScale(iii)=[];
                T(:,iii)=[];
                ii=ii+1;
                deleted = 20;
            end
        end
    end
    if maxChange<ThresholdEm
        if deleted > 0
            deleted = deleted -1;
        end
        if deleted==0
            break
        end
    end
allocationScore = sum(indices==(1:m));
allocationScore = mean(allocationScore(sum(indices==(1:m))>maxUsersPerCluster)-maxUsersPerCluster);
if isnan(allocationScore)
    allocationScore = 0;
end
reliabilityScore = sum(maxSINR > db2mag(0))/n + 1;
end