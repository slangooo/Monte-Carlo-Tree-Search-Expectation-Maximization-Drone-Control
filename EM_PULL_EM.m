function [allocationScore,reliabilityScore] = EM_PULL_EM (x,m,sigmaMax,maxUsersPerCluster,maxIter,SINR_THRESHOLD)
%Here the contribution of x to parameter update of drone j is determined by
ThresholdEm = 0.00001;
ThresholdDel = 0.5;
dSigmaMax=1;
sigmaMin = 1;
dMuMax = 1;
% sigmaMax = 1000;
AxisGranularity=4; % number of axis configurations to be tested
n=size(x,1); %number of test points
mu = x(chooseAxis(x,m,AxisGranularity),:);% this function return the indices of points to be used as initial means
initialSigma = 1;
sigma = [initialSigma initialSigma].*eye(2); %initialize sigma
sigma = repmat(sigma,1,1,m);
% softmaxScale = 40;
sigmaValPrev=initialSigma*ones(1,m);
SINR_desired_threshold = db2mag(20); %linear
% maxUsersPerCluster = 500;
deleted = 20;
iter = 0;
muNew = mu;
% maxIter = 500;
Pj = ones(1,m)/m; %initializing Pj
while iter < maxIter
    iter= iter + 1;
    input('')
    hold off;
    plotNetwork(x,sigma,mu,[min(x(:,1)) max(x(:,1))],[min(x(:,2)) max(x(:,2))]);
    muT=mu;
    sigmaT=sigma;
    [SINR, indices] = getSINR2(m,x,mu,sigma);
    usersPerCluster = sum(indices);
    TMask = logical( SINR < SINR_desired_threshold);
    for j = 1:m
        T(:,j)=Pj(j).*mvnpdf(x,mu(j,:),sigma(:,:,j));%obtain posterior probabilities (E)
    end
%     T = TMask.*T +  ~TMask.*indices.*T;
    
    T= T./sum(T,2); %normalize posterior probabilities
    Pj =sum(T,1);
    for j=1:m
        muNew(j,:)=[T(:,j)'*x(:,1) T(:,j)'*x(:,2)]./Pj(j);
        dMu = muNew(j,:) - mu(j,:);
        if norm(dMu) > dMuMax
            mu(j,:) = mu(j,:) + dMuMax * dMu./norm(dMu);
        else
            mu(j,:) = muNew(j,:);
        end
        Tx=T;
        sigmaVal= max(1,min(sigmaMax,(Tx(:,j)'*(x(:,1)-mu(j,1)).^2 + Tx(:,j)'*(x(:,2)-mu(j,2)).^2)/(Pj(j)))); %average of variances
        sigmaVal=min(sigmaVal, sigmaValPrev(j)+dSigmaMax);
        if ((sigmaVal < sigmaValPrev(j)) || (sigmaVal > sigmaValPrev(j) && usersPerCluster(j) < maxUsersPerCluster )) 
          
            sigma(:,:,j)= [sigmaVal 0; 0 sigmaVal];
            sigmaValPrev(j)=sigmaVal;
        elseif usersPerCluster(j) > maxUsersPerCluster && sigma(1,1,j) > sigmaMin && maxChange<ThresholdDel
            sigma(:,:,j)= sigma(:,:,j).*0.98;
        end
    end
    Pj = Pj./n;
    muT= abs(mu-muT);
    sigmaT=abs(sigma-sigmaT);
    maxChange=max(max(muT,[],'all') ,max(sigmaT,[],'all'));
    if (maxChange<ThresholdDel || iter == maxIter - 40)
        ii=0;
        for i=1:m
            iii=i-ii;
%             overlaps = circleOverlapArea(vecnorm(mu(iii,:)-mu,2,2),sqrt(sigma(1,1,iii)./1000),sqrt(squeeze(sigma(1,1,:))./1000));
%             overlaps(isnan(overlaps))=[];
                sqrt(squeeze(-sigma(1,1,iii) + sigma(1,1,:)))
                vecnorm(mu(iii,:)-mu,2,2)
            if any(vecnorm(mu(iii,:)-mu,2,2)<sqrt(squeeze(-sigma(1,1,iii) + sigma(1,1,:)))) || sigma(1,1,iii) < sigmaMin
                m=m-1;
                mu(iii,:)=[];
                muT(iii,:)=[];
                sigmaT(:,:,iii)=[];
                sigma(:,:,iii)=[];
                Pj(iii)=[];
%                 sigmaScale(iii)=[];
                T(:,iii)=[];
                sigmaValPrev(iii)=[];
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
reliabilityScore = sum(SINR > db2mag(SINR_THRESHOLD))/n;
end