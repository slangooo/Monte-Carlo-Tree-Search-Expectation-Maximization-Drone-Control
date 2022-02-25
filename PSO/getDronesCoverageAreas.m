function areas = getDronesCoverageAreas (mu,r,points, smallSteps, largeSteps)
m = size(mu,1);
areas = zeros(1,m);
for k =1:m
    area=0;
    for i=1:smallSteps-1
        for j=1:smallSteps-1
            if sqrt((mu(k,1)-points(i,j,1))^2+(mu(k,2)-points(i,j,2))^2)< r(k)
%             if norm(mu(k,:)- reshape(points(i,j,:),1,2)) < r(k)
                area = area +1;
            end
        end
    end
    areas(k) = area;
end