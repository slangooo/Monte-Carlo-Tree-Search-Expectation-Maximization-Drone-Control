function  [reliabilitiesPerIterationMcts, numberOfStatesVisitedMcts,energyEfficiencyPerIteration, numberOfServedUsersPerIteration, EmIter] = runMcts2 (x, parameters)
totalBandwidth =  parameters.totalBandwidth;
requiredRate =  parameters.requiredRate;
spectralEfficiency =  parameters.spectralEfficiency;
capacityBS = totalBandwidth.* spectralEfficiency;
numberOfUsersPerCluster = capacityBS./ requiredRate;
maxUsersPerCluster = numberOfUsersPerCluster;
bandwidthPerUser = totalBandwidth./maxUsersPerCluster;

m_initial = parameters.numberOfDrones;
% m_initial = ceil(size(x,1)/maxUsersPerCluster)+1;

global maxScore_K_05_N_5
global statesHistory_K_05_N_5
global childHistory_K_05_N_5
global reliabilities_K_05_N_5
global maxSigma
global maxAllocationScore_K_05_N_5;
global maxEfficiencyScore_K_05_N_5;

global EmIter
global EmIterCounter
EmIter = 0;
EmIterCounter = 0;

childHistory_K_05_N_5=[];
statesHistory_K_05_N_5 = zeros(1,m_initial);
maxScore_K_05_N_5 =0;
reliabilities_K_05_N_5 = [];
maxSigma =[];
maxAllocationScore_K_05_N_5 = 0;
maxEfficiencyScore_K_05_N_5 = 0;

softmaxScale = 4;
maxIter = parameters.maxIterMcts;
step = parameters.mctsStep;
m_initial = parameters.numberOfDrones;

desiredLinkReliability = parameters.linkReliabilityThreshold;
sinr_threshold =  parameters.sinr_threshold;
dronePower =  parameters.dronePower;
sigma = 1.*dronePower.*ones(1,m_initial);
n= size(x,1);
numberOfStatesVisitedMcts = zeros(maxIter,1);
sinrThresholds = -5:22;

reliabilitiesPerIterationMcts = zeros(length(sinrThresholds),maxIter);
energyEfficiencyPerIteration = zeros(maxIter,1);
numberOfServedUsersPerIteration = zeros(maxIter,1);

numberOfUsers = size(x,1);

root= mcts_node_rewards_F_noLoss(x,sigma,softmaxScale,maxUsersPerCluster,200,step,sinr_threshold,dronePower,numberOfUsers, desiredLinkReliability,0);

for iteration = 1:maxIter
    iterMcts = iteration
    root.update(x,root.numberOfServedUsers/root.numberOfUsers);
    numberOfStatesVisitedMcts(iteration) = EmIter;
    if ~isempty(reliabilities_K_05_N_5)
        reliabilitiesPerIterationMcts(:,iteration) =sum(10.*log10(reliabilities_K_05_N_5)>sinrThresholds)/n;
        energyEfficiencyPerIteration(iteration) = (bandwidthPerUser.*sum(log2(reliabilities_K_05_N_5 + 1)))./sum(maxSigma);
        numberOfServedUsersPerIteration(iteration) = maxAllocationScore_K_05_N_5;
    else
        reliabilitiesPerIterationMcts(:,iteration) =0;
    end
    
end
