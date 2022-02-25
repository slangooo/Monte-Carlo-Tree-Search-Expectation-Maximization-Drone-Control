function area = getOverlappingArea (points,mu,r)
area=0;
for i=1:size(points,1)
    for j=1:size(points,2)
        if sqrt((mu(1)-points(i,j,1))^2 + (mu(2)-points(i,j,2))^2)<r
%         if norm(mu- reshape(points(i,j,:),1,2)) < r
            area = area +1;
        end
    end
end