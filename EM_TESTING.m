x = [mvnrnd([-5 -5],reshape([80 50],1,2).*eye(2),100); mvnrnd([-5 9],reshape([20 35],1,2).*eye(2),100);mvnrnd([8 -2],reshape([35 20],1,2).*eye(2),100)];
x = [unifrnd(0,1000,400,2)];
x = [0 0; 50 50];
scatter(x(:,1),x(:,2))
m=4;
sigma_in=1.*100.*ones(1,m);
softmaxScale= 4;
maxUsersPerCluster=100;
maxIter=100;
sinr_threshold=2;
[ efficiencyScore,allocationScore,reliabilityScore, maxSINR, numberOfServedUsers,iter] = EM_PULL_CNST_SIGMA1...
                                                                (x,m,sigma_in,softmaxScale,maxUsersPerCluster,maxIter,sinr_threshold);