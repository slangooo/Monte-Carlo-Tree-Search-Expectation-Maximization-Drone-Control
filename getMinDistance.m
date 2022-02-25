function [minDistance] =  getMinDistance(maxMu)
minDistance = sqrt(sum((maxMu(1,:)-maxMu(2,:)).^2));
for ii =1:size(maxMu,1)
    for jj = 1:size(maxMu,1)
        if ii~=jj
            tmp = sqrt(sum((maxMu(ii,:)- maxMu(jj,:)).^2));
            if tmp < minDistance
                minDistance = tmp;
            end
        end
    end
end