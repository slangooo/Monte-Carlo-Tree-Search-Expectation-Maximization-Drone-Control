function SINR = getSINR(m,x,mu,sigma)
h=10; %height of drone
pathLossExponent=2;
refPower= -40; %reference power at 1m db
Noise = -95;% receiver noise power dBw;
n=size(x,1);
SNR_Thresh = 10; %SNR value to be ensured for all users
%%%
% % R=zeros(m,1);
% % K = zeros(n,m);
%get SEDs between all drones and users
% % for j=1:m
% %     K(:,j)=sqrt(sum((x-mu(j,:)).^2,2)+h^2).^pathLossExponent./sigma(1,1,j);
% % end
% % [V,Allocs]=min(K,[],2); %allocations of users to drones
% % binAllocsInv=logical(full(ind2vec(Allocs'))');
% % while size(binAllocsInv,2) < m
% %     binAllocsInv =logical([binAllocsInv zeros(n,1)]);
% % end

% % for j=1:m
% %     [V,I]=max(K(Allocs==j,j));
% %     V=x(Allocs==j,:);
% %     V=V(I,:);
% %     %get distance to fathest user
% %     R(j)= norm((V-mu(j,:)));
% % end
% After associating users with drone that has lowest SED, get drone
% transmit power based SNR to farthest user
% % TPower=getPower(R, pathLossExponent,refPower,Noise,h,SNR_Thresh);
%Or get TPower directly proportional to sigma
TPower= squeeze(deal(sigma(1,1,:)./1000));
TPower = sigma./1000;
indices = (1:m)'==(1:m);
%distances between each user and each drone
distancesN = sqrt((mu(:,1)'-x(:,1)).^2+ (mu(:,2)'-x(:,2)).^2 +h^2); 
distancesN=distancesN.^-pathLossExponent;
RPower=distancesN.*TPower.*db2mag(refPower);
SINR =[];
for i = 1: m
  	SINR_T = RPower(repmat(indices(i,:),n,1))./(db2mag(Noise)+sum(reshape(RPower(~repmat(indices(i,:),n,1)),n,m-1),2));
    SINR = cat(2,SINR,SINR_T);
end
% RTotal=sum(log2(1+SINR)); %total SINR
% PTotal= sum(TPower); %total Power
% EE=RTotal/PTotal;
end
