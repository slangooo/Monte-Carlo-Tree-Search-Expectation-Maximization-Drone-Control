function U1 = getFirstUtility(ro,density1, density2, maxUsersPerCluster, xBoundary ,smallSteps, largeSteps,m)
segsPerArea = (smallSteps -1)/(largeSteps-1);
U1 = 0;

for i=0:largeSteps-2
    for j=0:largeSteps-2
        UU =0;
        if i*segsPerArea>=xBoundary(2)/xBoundary(3)*(smallSteps-1)
            density = density2;
        else
            density = density1;
        end
        for k=1:m
            UU = UU + ro(i+1,j+1,k)*maxUsersPerCluster ;%-0.1;
        end
        UU = UU - density;
        if UU > 0
            UU = 0;
        end
        U1 = U1 + UU;
    end
end