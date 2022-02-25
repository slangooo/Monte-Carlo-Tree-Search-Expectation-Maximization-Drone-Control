function [mu,sigma,m,Pj,EE,Iscore,TurnedOff]= EM_Monte(x, inSigma, mu, Pj, m, scaleSigma)
MaxSigma=10;
MinSigma=2;
AxisGranularity=4;
ThresholdEm=0.00001;
dSigmaMax=0.2;
pathLossExponent=2;
FLAG=0;
maxChange=1;
n=size(x,1); %number of test points
refPower = -40; %reference power at 1m db
Noise= -95; % receiver noise power dBw;
h=10; %height m
maxUsersPerCluster=100; %number of maximum users per cluster
K=zeros(n,m);
added=10;
global depthFactor;
depthFactor=depthFactor;
% depthFactor=0.9; %Depth Vs. Breadth
TurnedOff=[];
mInitial=m;
%TODOOOOO preallocation instead of dynamic!
if isempty(mu)
    mu = x(chooseAxis(x,m,AxisGranularity),:);% this function return the indices of points to be used as initial means
end
sigma = inSigma; %initialize sigma
if isempty(Pj)
    Pj = ones(1,m)/m; %initializing Pj
end
T=zeros(n,m); %=Pjt
sigmaValPrev=zeros(1,m);
for j=1:m
    sigmaValPrev(j)=inSigma(1,1,j);
end
if ~isempty(scaleSigma) %1-scale 0-no constraint -1 - dont change
    for j=1:m
        if scaleSigma(j)==1
            sigma(:,:,j)=sigma(:,:,j).*0.9;
        end
    end
end
while 1
    muT=mu;
    sigmaT=sigma;
    for j = 1:m
        T(:,j)=Pj(j).*mvnpdf2(x,mu(j,:),sigma(:,:,j));%obtain posterior probabilities (E)
    end
    T= T./sum(T,2); %normalize posterior probabilities
    Pj =sum(T,1);
    for j=1:m
        if FLAG==0 || FLAG==1
            mu(j,:)=[T(:,j)'*x(:,1) T(:,j)'*x(:,2)]./Pj(j);
        end
        if FLAG==0 || FLAG==2
            if  j>size(scaleSigma,2) || scaleSigma(j)==0
                sigmaVal= max(MinSigma,min(MaxSigma,(T(:,j)'*(x(:,1)-mu(j,1)).^2 + T(:,j)'*(x(:,2)-mu(j,2)).^2)/(2*Pj(j)))); %average of variances
                sigmaVal=min(sigmaVal, sigmaValPrev(j)+dSigmaMax);
                sigma(:,:,j)= [sigmaVal 0; 0 sigmaVal];
                sigmaValPrev(j)=sigmaVal;
            end
        end
    end
    
    
    for j=1:m
        K(:,j)=sqrt(sum((x-mu(j,:)).^2,2)+h^2).^pathLossExponent./sigma(1,1,j);
    end
    [V,Allocs]=min(K,[],2); %allocations of users to drones
    binAllocsInv=logical(full(ind2vec(Allocs'))');
    while size(binAllocsInv,2) < m
        binAllocsInv =logical([binAllocsInv zeros(n,1)]);
    end
    usersDoNotExceed=1;
    for j=1:m
        if sum(binAllocsInv(:,j))>maxUsersPerCluster
            usersDoNotExceed=0;
            added=added-1;
            if maxChange<ThresholdEm*1e1 && added <=0
                added=10;
                m=m+1;
                mu=[mu;mu(j,:)+0.5];
                muT=[muT;mu(j,:)+0.5];
                mu(j,:)= mu(j,:)-0.5;
                muT(j,:)= muT(j,:)-0.5;
                sigma=cat(3,sigma,sigma(:,:,j));
                sigmaT=cat(3,sigmaT,sigma(:,:,j));
                for iii=1:m
                    val=max(sigma(1,1,iii)/2,MinSigma);
                    sigma(:,:,iii)=[val val].*eye(2);
                    sigmaValPrev(iii)=sigma(1,1,iii);
                end
                
                %                 sigma=cat(3,sigma,[2 2].*eye(2));
                %                 sigmaT=cat(3,sigmaT,[2 2].*eye(2));
                K=[K zeros(n,1)];
                Pj(j)=Pj(j)/2;
                Pj=[Pj Pj(j)];
                T=[T zeros(n,1)];
                sigmaValPrev=[sigmaValPrev sigma(1,1,j)];
                sigmaValPrev(j)=sigma(1,1,j);
            end
        end
    end
    
    muT= abs(mu-muT);
    sigmaT=abs(sigma-sigmaT);
    maxChange=max(max(muT,[],'all') ,max(sigmaT,[],'all'));
    Pj = Pj./n;
    if maxChange<ThresholdEm && usersDoNotExceed
        break
    end
end
jj=0;
for j=1:m
    jjj=j-jj;
    if sum(binAllocsInv(:,jjj)) < 1 %%%maybe remove
        m=m-1;
        mu(jjj,:)=[];
        muT(jjj,:)=[];
        sigmaT(:,:,jjj)=[];
        sigma(:,:,jjj)=[];
        K(:,jjj)=[];
        Pj(jjj)=[];
        T(:,jjj)=[];
        sigmaValPrev(jjj)=[];
        jj=jj+1;
        binAllocsInv(:,jjj)=[];
        if j<=mInitial
            TurnedOff=[TurnedOff j];
        end
    end
end
R=zeros(m,1);
for j=1:m
    K(:,j)=sqrt(sum((x-mu(j,:)).^2,2)+h^2).^pathLossExponent./sigma(1,1,j);
end
[V,Allocs]=min(K,[],2); %allocations of users to drones
binAllocsInv=logical(full(ind2vec(Allocs'))');
for j=1:m
    [V,I]=max(K(Allocs==j,j));
    V=x(Allocs==j,:);
    V=V(I,:);
    R(j)= norm((V-mu(j,:)));
end

TPower=getPower(R, pathLossExponent,refPower,Noise,h); %based on minimum SNR
binAllocs =~binAllocsInv; % return one-hot reversed encoding of Allocs
distancesN = sqrt((mu(:,1)'-x(:,1)).^2+ (mu(:,2)'-x(:,2)).^2 +h^2); %distances between each user and each drone
distancesN=distancesN.^-pathLossExponent;
RPower=distancesN.*TPower'.*db2mag(refPower)./db2mag(Noise);
SINR=RPower(binAllocsInv')./sum(reshape(RPower(binAllocs),m-1,n)',2); %SINR per user
RTotal=sum(log2(1+SINR)); %total SINR
PTotal= sum(TPower); %total Power
EE=RTotal/PTotal;
Iscore=sum(binAllocs.*RPower);
Iscore= softmax((((Iscore-max(Iscore))/(depthFactor* min(Iscore)))'));
