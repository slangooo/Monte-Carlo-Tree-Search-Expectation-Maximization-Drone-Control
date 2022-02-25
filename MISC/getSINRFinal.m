function SINR = getSINRFinal (x,mu,sigma,h)
a = 9.61;
b = 0.16;
carrFreq = 2e9;
avgLossLOS = 1;
avgLossNLOS = 20;
totalBandwidth = 20e6;
requiredRate = 1e6;
spectralEfficiency = 2;
noisePower = -140;


capacityBS = totalBandwidth.* spectralEfficiency;
numberOfUsersPerCluster = capacityBS./ requiredRate;
userBandwidth = totalBandwidth / numberOfUsersPerCluster;


distances2D = sqrt((mu(:,1)'-x(:,1)).^2+ (mu(:,2)'-x(:,2)).^2);
distances3D = sqrt(distances2D.^2 + h.^2);
PLOS = getProbLOS(h,distances2D,a,b);
PL = getPathLoss (carrFreq, distances3D, PLOS, avgLossLOS, avgLossNLOS);
Pr = getReceivedPower(PL, sigma, userBandwidth, totalBandwidth);

% SNR = getSNR (Pr, noisePower);

sigma(sigma==0) = [];
m=sum(sigma~=0);
indices = (1:m)'==(1:m);
n = size(x,1);
SINR = zeros(n,m);
for i = 1: m
    %   	SINR_T = Pr(repmat(indices(i,:),n,1))./(10.^(noisePower/10)+sum(reshape(Pr(~repmat(indices(i,:),n,1)),n,m-1),2));
    SINR_T = Pr(:,i);%./(10.^(noisePower/10)+sum(reshape(Pr(~repmat(indices(i,:),n,1)),n,m-1),2));
    SINR_D = zeros(n,1);
    for j=1:m
        if j == i
            continue
        end
        SINR_D = SINR_D + Pr(:,j);
    end
    SINR(:,i)= SINR_T./(10.^(noisePower/10)+SINR_D);
%     SINR(:,i)= SINR_T;
end