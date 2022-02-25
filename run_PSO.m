a = 9.61;
b = 0.16;
carrFreq = 2e9;
avgLossLOS = 1;
avgLossNLOS = 20;
totalBandwidth = 20e6;
requiredRate = 1e6;
spectralEfficiency = 2;
noisePower = -140;
capacityBS = totalBandwidth.* spectralEfficiency;
numberOfUsersPerCluster = capacityBS./ requiredRate;
userBandwidth = totalBandwidth / numberOfUsersPerCluster;

maxIterations = 200;

sinr_threshold = 2;
SNRThresh = 15;
linkReliabilityThreshold = 0.95;

dronePower = 500;
numberOfDrones = 20;

xBoundary =[0 200 1000];
yBoundary =[0 1000];
largeSteps = 6;
smallSteps = 51;

numberOfUsers = [200 300];

c1 = 2; %local
c2 = 1; %global
phi = 0.1;
phi1Range = [0 1]; %local
phi2Range = [0 1]; %global

sinrThresholds = -5:22;%%%%%%%%%%



reliabilitiesPerIteration = zeros(length(sinrThresholds),maxIterations);

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
    mus(i,:,:)=unifrnd(0,xBoundary(3),m_initial,2);
    heights(i,:) = unifrnd(10,200,1,m_initial);
    [U1,U2,U3,maxSINRParticle,associations,reliability] = getUtilitiesFromParticle (x, reshape(mus(i,:,:),m_initial,2), heights(i,:), sigma, maxUsersPerCluster,...
        sinr_threshold, SNRThresh, noisePower, userBandwidth, totalBandwidth, avgLossNLOS,...
        avgLossLOS, carrFreq,a,b,spectralEfficiency, density1, density2, xBoundary, yBoundary,...
        smallSteps, largeSteps);
    Ulocal(i) = U1;
    Wlocal(i,:,:) =  [reshape(mus(i,:,:),m_initial,2) reshape(heights(i,:),m_initial,1)];
    if i==1
        Uglobal = U1;
        Wglobal = [reshape(mus(i,:,:),m_initial,2) reshape(heights(i,:),m_initial,1)];
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

phi1 = unifrnd(phi1Range(1),phi1Range(2));
phi2 = unifrnd(phi2Range(1),phi2Range(2));
Vprev=zeros(numberOfParticles,m_initial,3);
% for i=1:numberOfParticles
%     V = phi.*squeeze(Vprev(i,:,:)) + c1*phi1.*squeeze(Wlocal(i,:,:)-W(i,:,:))+ c2*phi2.*(Wglobal-squeeze(W(i,:,:)));
%     Vprev(i,:,:) = V;
%     W(i,:,:) = squeeze(W(i,:,:)) +V;
% end
maxAssoc=[];
Uindex =1;
U=zeros(1,3);
while 1
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
            smallSteps, largeSteps);
        if Ulocal(i)> U(Uindex)
            if (Uindex~=1 && U(1) <=0) || Uindex ==1
                Ulocal(i) = U(Uindex);
                Wlocal(i,:,:) =  W(i,:,:);
            end
        end
        if Uglobal > Ulocal(i)
            if (Uindex~=1 && U(1) <=0) || Uindex ==1
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
%             Uglobal = 
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
    if iterations > maxIterations
        break
    end
    reliabilitiesPerIteration(:,iterations)=sum(10*log10(maxSINR)>sinrThresholds)/n;
end

% figure
% hold on
% for j = 1:m_initial
%     a=20; % horizontal radius
%     b=a; % vertical radius
%     x0=Wglobal(j,1); % x0,y0 ellipse centre coordinates
%     y0=Wglobal(j,2);
%     t=-pi:0.01:pi;
%     xp=x0+a*cos(t);
%     yp=y0+b*sin(t);
%     plot(xp,yp)
% end
% sum(maxAssoc==(1:m_initial))
% sum(10*log10(maxSINR)>sinr_threshold)/n
% sum(10*log10(maxSINR)>sinrThresholds)/n
% sinrThresholds
