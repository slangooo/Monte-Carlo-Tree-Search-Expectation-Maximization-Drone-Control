function [reliableUsersPerIteration,reliabilitiesPerIteration,x, energyEfficiencyPerIteration, numberOfServedUsersPerIteration] = runPso(parameters)

a =parameters.a;
b = parameters.b;
carrFreq = parameters.carrFreq;
avgLossLOS = parameters.avgLossLOS;
avgLossNLOS =  parameters.avgLossNLOS;
totalBandwidth =  parameters.totalBandwidth;
requiredRate =  parameters.requiredRate;
spectralEfficiency =  parameters.spectralEfficiency;
noisePower =  parameters.noisePower;

capacityBS = totalBandwidth.* spectralEfficiency;
numberOfUsersPerCluster = capacityBS./ requiredRate;
userBandwidth = totalBandwidth / numberOfUsersPerCluster;

maxIterations =  parameters.maxIterations;

sinr_threshold =  parameters.sinr_threshold;
SNRThresh =  parameters.SNRThresh;
linkReliabilityThreshold =  parameters.linkReliabilityThreshold;

dronePower =  parameters.dronePower;
numberOfDrones =  parameters.numberOfDrones;
% numberOfDrones = ceil(sum(parameters.numberOfUsers)/numberOfUsersPerCluster)+1;


xBoundary = parameters.xBoundary;
yBoundary = parameters.yBoundary;
largeSteps =  parameters.largeSteps;
smallSteps =  parameters.smallSteps;

numberOfUsers =  parameters.numberOfUsers;

c1 =  parameters.c1; %local
c2 =  parameters.c2; %global
% phi =  parameters.phi;
phi1Range =  parameters.phi1Range; %local
phi2Range =  parameters.phi2Range; %global

sinrThresholds = -5:22;%%%%%%%%%%



reliabilitiesPerIteration = zeros(length(sinrThresholds),maxIterations);
energyEfficiencyPerIteration = zeros(maxIterations,1);
reliableUsersPerIteration= zeros(maxIterations,1);
numberOfServedUsersPerIteration = zeros(maxIterations,1);

m_initial = numberOfDrones;
sigma = 1.*dronePower.*ones(1,m_initial);

y = unifrnd(0,yBoundary(2),sum(numberOfUsers),1);
x=unifrnd(0,xBoundary(2),numberOfUsers(1),1);
x = [x ; unifrnd(xBoundary(2),xBoundary(3),numberOfUsers(2),1)];
x = [x y];

maxUsersPerCluster = numberOfUsersPerCluster;

numberOfLargeXSegments = [xBoundary(2)/xBoundary(3).*(largeSteps-1) (xBoundary(3)- xBoundary(2))/xBoundary(3).*(largeSteps-1)];
% density = (numberOfUsers(1)/(largeSteps-1)^2);
density1 = numberOfUsers(1)/(numberOfLargeXSegments(1)*(largeSteps-1));
density2 = numberOfUsers(2)/(numberOfLargeXSegments(2)*(largeSteps-1));

[subX, subY, X, Y] = createMesh (xBoundary, yBoundary, smallSteps, largeSteps);
gridPoints = cat(3,subX,subY);


numberOfParticles = 3 * m_initial;

n= size(x,1);

Ulocal = zeros(numberOfParticles,1);
Uglobal = 0;
Wlocal = zeros(numberOfParticles,m_initial,3);
W = zeros(numberOfParticles,m_initial,3);

iterations = 0;
mus = zeros(numberOfParticles,m_initial,2);
heights = zeros(numberOfParticles,m_initial);
for i=1:numberOfParticles
    mus(i,:,:)=unifrnd(-50,xBoundary(3)+50,m_initial,2);
    heights(i,:) = unifrnd(60,90,1,m_initial);
    [U1,U2,U3,maxSINRParticle,associations,reliability] = getUtilitiesFromParticle (x, reshape(mus(i,:,:),m_initial,2), heights(i,:), sigma, maxUsersPerCluster,...
        sinr_threshold, SNRThresh, noisePower, userBandwidth, totalBandwidth, avgLossNLOS,...
        avgLossLOS, carrFreq,a,b,spectralEfficiency, density1, density2, xBoundary, yBoundary,...
        smallSteps, largeSteps, gridPoints);
    Ulocal(i) = U1;
    Wlocal(i,:,:) =  [reshape(mus(i,:,:),m_initial,2) reshape(heights(i,:),m_initial,1)];
    if i==1
        Uglobal = U1;
        Wglobal = [reshape(mus(i,:,:),m_initial,2) reshape(heights(i,:),m_initial,1)];
        maxSINR = maxSINRParticle;
        maxAssoc = associations;
    end
    if Uglobal > U1
        Uglobal = U1;
        Wglobal = [reshape(mus(i,:,:),m_initial,2) reshape(heights(i,:),m_initial,1)];
        maxSINR = maxSINRParticle;
        maxAssoc = associations;
        maxReliab = reliability;
    end
end

W = Wlocal;

% muPrev = mus;
% heightsPrev = heights;
% Wprev = cat(3,muPrev,heightsPrev);

% phi1 = unifrnd(phi1Range(1),phi1Range(2));
% phi2 = unifrnd(phi2Range(1),phi2Range(2));
Vprev=zeros(numberOfParticles,m_initial,3);
% for i=1:numberOfParticles
%     V = phi.*squeeze(Vprev(i,:,:)) + c1*phi1.*squeeze(Wlocal(i,:,:)-W(i,:,:))+ c2*phi2.*(Wglobal-squeeze(W(i,:,:)));
%     Vprev(i,:,:) = V;
%     W(i,:,:) = squeeze(W(i,:,:)) +V;
% end
Uindex =1;
U=zeros(1,3);

%%%%%plotting
% colors = ['c','y','m','r','g','b','w','k'];
% axesH = axes('XLim',[0 1000], 'YLim', [0 1000], 'ZLim', [0 200], ...
%     'NextPlot', 'add');
% grid on;
% view(3);
% plot3(axesH,x(:,1),x(:,2),zeros(length(x),2),'.');
% dataPlot = plot3(axesH,W(i,:,1),W(i,:,2),W(i,:,3),strcat('o',colors(mod(i,8)+1)));
% figure
% plot3(Wglobal(:,1),Wglobal(:,2),Wglobal(:,3),'o')
%%%%%%%%

while 1
    phi = parameters.phiMax - (parameters.phiMax - parameters.phiMin)*iterations/parameters.maxIterations;
    %%%%%
%     pause(0.05)
%     delete(dataPlot);
% %     dataPlot = plot3(Wglobal(:,1),Wglobal(:,2),Wglobal(:,3),'o');
%     for i=1:numberOfParticles
%         dataPlot(i) = plot3(axesH,W(i,:,1),W(i,:,2),W(i,:,3),strcat('o',colors(mod(i,8)+1)));
%     end
%     drawnow;
    %%%%
    iterPso = iterations
    iterations = iterations + 1;
    for i=1:numberOfParticles
        phi1 = unifrnd(phi1Range(1),phi1Range(2));
        phi2 = unifrnd(phi2Range(1),phi2Range(2));
        V = phi.*squeeze(Vprev(i,:,:)) + c1*phi1.*squeeze(Wlocal(i,:,:)-W(i,:,:))+ c2*phi2.*(Wglobal-squeeze(W(i,:,:)));
        Vprev(i,:,:) = V;
        W(i,:,:) = squeeze(W(i,:,:)) +V;
    end
    for i=1:numberOfParticles
        mus(i,:,:)=squeeze(W(i,:,1:2));
        heights(i,:) = squeeze(W(i,:,3));
        [U(1),U(2),U(3),maxSINRParticle,associations,reliability] = getUtilitiesFromParticle (x, reshape(mus(i,:,:),m_initial,2), heights(i,:), sigma, maxUsersPerCluster,...
            sinr_threshold, SNRThresh, noisePower, userBandwidth, totalBandwidth, avgLossNLOS,...
            avgLossLOS, carrFreq,a,b,spectralEfficiency, density1, density2, xBoundary, yBoundary,...
            smallSteps, largeSteps, gridPoints);
        if Ulocal(i)> U(Uindex)
            if (Uindex==3 && U(1) <=0 && U(2) <=-linkReliabilityThreshold*n)...
                    ||Uindex ==2 && U(1) <=0 ||Uindex ==1
                Ulocal(i) = U(Uindex);
                Wlocal(i,:,:) =  W(i,:,:);
            end
        end
        if Uglobal > Ulocal(i)
            if (Uindex==3 && U(1) <=0 && U(2) <=-linkReliabilityThreshold*n)...
                    ||Uindex ==2 && U(1) <=0 ||Uindex ==1
                Uglobal = Ulocal(i);
                Wglobal = squeeze(Wlocal(i,:,:));
                maxSINR = maxSINRParticle;
                maxAssoc = associations;
                maxReliab = reliability;
            end
        end
    end
    if Uglobal <=0
%         if Uindex ==1
%             Uglobal = 0;
%         end
        Uindex = 2;
    end
    if Uglobal <=-linkReliabilityThreshold*n
        Uindex=3;
    end
    Wglobal(abs(Wglobal(:,3))<10,3) = 10;
    Wglobal(Wglobal(:,3)<-10,3) = -Wglobal(Wglobal(:,3)<-10,3);
    if  Uglobal <= -n
        break
    end
    Uglobal
    if iterations > maxIterations
        break
    end
    reliabilitiesPerIteration(:,iterations:end)=repmat(sum(10*log10(maxSINR)>sinrThresholds)'/n,1,maxIterations-iterations+1);
    energyEfficiencyPerIteration(iterations:end) = (userBandwidth.*sum(log2(maxSINR + 1)))./sum(sigma);
    usersExceeding = sum(maxAssoc == (1:numberOfDrones)) - maxUsersPerCluster;
    usersExceeding(usersExceeding<0)=0;
    numberOfServedUsersPerIteration(iterations:end) = n-sum(usersExceeding,2);
    
    usersExceeding = sum(((maxSINR > 10.^(sinr_threshold./10)).*maxAssoc)==(1:numberOfDrones))...
    - maxUsersPerCluster;
    usersExceeding(usersExceeding<0)=0;
    reliableUsers = (sum(maxSINR > 10.^(sinr_threshold./10))- sum(usersExceeding));
    reliableUsersPerIteration(iterations:end) = reliableUsers;
end
