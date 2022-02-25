function ro = getDronesSubAreasOverlap(mu ,points, largeSteps, smallSteps, r, areas)
m=size(mu,1);
ro = zeros(largeSteps-1,largeSteps-1,m);
segsPerArea = (smallSteps -1)/(largeSteps-1);
for k=1:m
    if areas(k)==0
        continue
    end
    for i=0:largeSteps-2
        for j=0:largeSteps-2
            subAreaPoints = points(i*segsPerArea+1:(i+1)*segsPerArea,j*segsPerArea+1:(j+1)*segsPerArea,:);
            subAreaOverlap=getOverlappingArea(subAreaPoints,mu(k,:),r(k));
            ro(i+1,j+1,k)=subAreaOverlap/areas(k);
        end
    end
end