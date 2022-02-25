a = 9.61;
b = 0.16;
carrFreq = 2e9; 
avgLossLOS = 1;
avgLossNLOS = 20;
totalBandwidth = 20e6;
requiredRate = 1e6;
spectralEfficiency = 2;
noisePower = -120;
capacityBS = totalBandwidth.* spectralEfficiency;
numberOfUsersPerCluster = capacityBS./ requiredRate;
userBandwidth = totalBandwidth / numberOfUsersPerCluster;

x=unifrnd(0,500,500,2);
maxUsersPerCluster = numberOfUsersPerCluster;
sinr_threshold = -2;
SNRThresh = 30;
linkReliabilityThreshold = 0.95;
mu = [0 10; 20 100; 20 40];

h = [10 50 40];
sigma = 5000.*ones(1,m);

density1 = 2e-3;
density2 = 2e-3;
xBoundary =[0 300 500];
yBoundary =[0 500];

largeSteps = 11;
smallSteps = 101;

[U1,U2,U3] = getUtilitiesFromParticle (x, mu, h, sigma, maxUsersPerCluster, sinr_threshold, SNRThresh, noisePower, userBandwidth, totalBandwidth, avgLossNLOS, avgLossLOS, carrFreq,a,b,spectralEfficiency, density1, density2, xBoundary, yBoundary, smallSteps, largeSteps);
% [subX, subY, X, Y] = createMesh (xBoundary, yBoundary, smallSteps, largeSteps);
% 
% r = getCoverageRadius(sigma./1000,SNRThresh,noisePower, h, userBandwidth, totalBandwidth, avgLossNLOS, avgLossLOS, carrFreq,a,b); 
% gridPoints = cat(3,subX,subY);
% droneCoverageAreas = getDronesCoverageAreas (mu,r,gridPoints, smallSteps, largeSteps);
% ro = getDronesSubAreasOverlap(mu ,gridPoints, largeSteps, smallSteps, r, droneCoverageAreas);
% U1 = getFirstUtility(ro,density1, density2, maxUsersPerCluster, xBoundary ,smallSteps, largeSteps,m);
% 
% 
% SINR = getSINRFinal (x,mu,sigma./1000,h);
% [maxSINR, indices] = max(SINR,[],2);
% U2 = -sum(maxSINR > 10.^(sinr_threshold./10));
% 
% userSpectralEfficiency = log2(1+maxSINR);
% averageSpectralEfficiency = mean(1./userSpectralEfficiency);
% U3 = -n + spectralEfficiency - 1/averageSpectralEfficiency;
% 
% 
figure(1)
plot(X,Y,'k')
hold on
plot(Y,X,'k')
% hold off
plot(subX,subY,':y')
% hold on
plot(subY,subX,':y')

for j = 1:m
    a=r(j); % horizontal radius
    b=a; % vertical radius
    x0=mu(j,1); % x0,y0 ellipse centre coordinates
    y0=mu(j,2);
    t=-pi:0.01:pi;
    xp=x0+a*cos(t);
    yp=y0+b*sin(t);
    plot(xp,yp)
end
% 
% 
% 
% 

%
% area1 = boundary(1,1) * boundary(1,2);
% area2 = (boundary(2,1) -boundary(1,1)) * boundary(2,2);
% numberOfAxisSegments = floor(sqrt(m));
%
% x = linspace(0, 100);
%
%
% userDensity1 = 2*10^-3;
% userDensity2 =