function [U1,U2,U3,maxSINR,associations,reliability] = getUtilitiesFromParticle (x, mu, h, sigma, maxUsersPerCluster, sinr_threshold, SNRThresh, noisePower, userBandwidth, totalBandwidth, avgLossNLOS, avgLossLOS, carrFreq,a,b,spectralEfficiency, density1, density2, xBoundary, yBoundary, smallSteps, largeSteps,gridPoints)
m = size(mu,1);
n=size(x,1);

% r = getCoverageRadius2(sigma./1000,SNRThresh,noisePower, h, userBandwidth, totalBandwidth, avgLossNLOS, avgLossLOS, carrFreq,a,b); 
% droneCoverageAreas = getDronesCoverageAreas (mu,r,gridPoints, smallSteps, largeSteps);
% ro = getDronesSubAreasOverlap(mu ,gridPoints, largeSteps, smallSteps, r, droneCoverageAreas);
% U1 = -1*getFirstUtility(ro,density1, density2, maxUsersPerCluster, xBoundary ,smallSteps, largeSteps,m);
% U1 = max(-1,U1);
% if U1<20
%     U1 = -1;
% end
U1 = 0; %Remove

SINR = getSINRFinal (x,mu,sigma./1000,h);
[maxSINR, indices] = max(SINR,[],2);
% [maxSNR, indices] = max(SNR,[],2);
associations = indices;
% usersExceeding = sum(associations == (1:m)) - maxUsersPerCluster;
usersExceeding = sum(((maxSINR > 10.^(sinr_threshold./10)).*associations)==(1:m))...
    - maxUsersPerCluster;
usersExceeding(usersExceeding<0)=0;

U2 = -(sum(maxSINR > 10.^(sinr_threshold./10))- sum(usersExceeding));
reliability = (sum(maxSINR > 10.^(sinr_threshold./10))- sum(usersExceeding))/n;
% U2 =  -sum(maxSINR > 10.^(sinr_threshold./10));

userSpectralEfficiency = log2(1+maxSINR);
averageSpectralEfficiency = mean(1./userSpectralEfficiency);
U3 = -n + spectralEfficiency - 1/averageSpectralEfficiency;

% figure(1)
% plot(X,Y,'k')
% hold on
% plot(Y,X,'k')
% % hold off
% plot(subX,subY,':y')
% % hold on
% plot(subY,subX,':y')
% 
% for j = 1:m
%     a=r(j); % horizontal radius
%     b=a; % vertical radius
%     x0=mu(j,1); % x0,y0 ellipse centre coordinates
%     y0=mu(j,2);
%     t=-pi:0.01:pi;
%     xp=x0+a*cos(t);
%     yp=y0+b*sin(t);
%     plot(xp,yp)
% end
