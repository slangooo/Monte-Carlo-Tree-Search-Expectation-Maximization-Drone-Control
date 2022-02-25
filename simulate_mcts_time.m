function [numberOfServedUsersPerIterationRandom, energyEfficiencyPerIterationRandom, numberOfStatesVisitedMcts, reliabilitiesPerIterationMcts, energyEfficiencyPerIteration, spectralEfficiencyPerIteration, numberOfServedUsersPerIteration, averageDronePower, minDistanceBetweenDronesPerIteration] =simulate_mcts_time(s_mobility,s_input,time_step,s_mobility_group,parameters)
v_t = 0:time_step:s_input.SIMULATION_TIME;
store= struct('v_x', zeros(1,1+s_input.SIMULATION_TIME/time_step),...
    'v_y',zeros(1,1+s_input.SIMULATION_TIME/time_step));
leaders=repmat(store,1,s_input.NB_NODES);
if ~isempty(s_mobility_group)
    cellGroup=struct2cell(s_mobility_group);
    Nmembers=cell2mat(cellGroup(1,:));
    cellGroup=[];
    members=repmat(store,1,sum(Nmembers));
end
countMembers=0;
for leaderIndex = 1:s_input.NB_NODES
    while s_mobility.VS_NODE(leaderIndex).V_TIME(end) - s_mobility.VS_NODE(leaderIndex).V_TIME(end-1)<0.1
        s_mobility.VS_NODE(leaderIndex).V_TIME(end)=[];
        s_mobility.VS_NODE(leaderIndex).V_POSITION_X(end)=[];
        s_mobility.VS_NODE(leaderIndex).V_POSITION_Y(end)=[];
    end
    leaders(leaderIndex).v_x = interp1(s_mobility.VS_NODE(leaderIndex).V_TIME,s_mobility.VS_NODE(leaderIndex).V_POSITION_X,v_t);
    leaders(leaderIndex).v_y = interp1(s_mobility.VS_NODE(leaderIndex).V_TIME,s_mobility.VS_NODE(leaderIndex).V_POSITION_Y,v_t);
    if isnan(leaders(leaderIndex).v_y(end))
        leaders(leaderIndex).v_x(end) = leaders(leaderIndex).v_x(end-1);
        leaders(leaderIndex).v_y(end) = leaders(leaderIndex).v_y(end-1);
    end
    if ~isempty(s_mobility_group)
        for groupIndex =1:s_mobility_group(leaderIndex).NB_NODES_group
            %         if numel(s_mobility.VS_NODE(leaderIndex).V_TIME)==2
            %             s_mobility_group(leaderIndex).VS_NODE(groupIndex).V_TIME(end)=[];
            %             s_mobility_group(leaderIndex).VS_NODE(groupIndex).V_POSITION_X(end)=[];
            %             s_mobility_group(leaderIndex).VS_NODE(groupIndex).V_POSITION_Y(end)=[];
            %         end
            members(groupIndex+countMembers).v_x=interp1(s_mobility_group(leaderIndex).VS_NODE(groupIndex).V_TIME,s_mobility_group(leaderIndex).VS_NODE(groupIndex).V_POSITION_X,v_t);
            members(groupIndex+countMembers).v_y=interp1(s_mobility_group(leaderIndex).VS_NODE(groupIndex).V_TIME,s_mobility_group(leaderIndex).VS_NODE(groupIndex).V_POSITION_Y,v_t);
            %         interp1(s_mobility_group(leaderIndex).VS_NODE(groupIndex).V_TIME,s_mobility_group(leaderIndex).VS_NODE(groupIndex).V_POSITION_Y,v_t)
        end
        countMembers=countMembers+s_mobility_group(leaderIndex).NB_NODES_group;
    end
end
leadersM=cell2mat(struct2cell(leaders));
if ~isempty(s_mobility_group)
    membersM=cell2mat(struct2cell(members));
end
vars = {'cellGroup', 'store', 'countMembers','members','leaders'};
clear(vars{:});
if ~isempty(s_mobility_group)
    totalMembers=sum(Nmembers);
end
timeIndexInitital = 1;
if ~isempty(s_mobility_group)
    x=[reshape(leadersM(:,timeIndexInitital,:),2,s_input.NB_NODES) reshape(membersM(:,timeIndexInitital,:),2,totalMembers)]';
else
    x=[reshape(leadersM(:,timeIndexInitital,:),2,s_input.NB_NODES)]';
end
x(all(x==0,2),:)=[];

timeSteps = length(v_t);

%%%%%%%%%%%%%%%MCTS
totalBandwidth =  parameters.totalBandwidth;
requiredRate =  parameters.requiredRate;
spectralEfficiency =  parameters.spectralEfficiency;
capacityBS = totalBandwidth.* spectralEfficiency;
numberOfUsersPerCluster = capacityBS./ requiredRate;
maxUsersPerCluster = numberOfUsersPerCluster;
bandwidthPerUser = totalBandwidth./maxUsersPerCluster;

m_initial = parameters.numberOfDrones;
global maxScore_K_05_N_5
global statesHistory_K_05_N_5
global childHistory_K_05_N_5
global reliabilities_K_05_N_5
global maxSigma
global maxAllocationScore_K_05_N_5;
global maxEfficiencyScore_K_05_N_5;
global EmIter
global EmIterCounter
global maxMu;
global firstVisitsCounter;

firstVisitsCounter=0;
maxMu = zeros(1,m_initial);
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
% timeSteps = parameters.timeStepsMcts;
step = parameters.mctsStep;
m_initial = parameters.numberOfDrones;

desiredLinkReliability = parameters.linkReliabilityThreshold;
sinr_threshold =  parameters.sinr_threshold;
dronePower =  parameters.dronePower;
sigma = 1.*dronePower.*ones(1,m_initial);
n= size(x,1);
numberOfStatesVisitedMcts = zeros(timeSteps,1);
sinrThresholds = -5:22;
averageDronePower = zeros(timeSteps,1);

reliabilitiesPerIterationMcts = zeros(length(sinrThresholds),timeSteps);
energyEfficiencyPerIteration = zeros(timeSteps,1);
spectralEfficiencyPerIteration = zeros(timeSteps,1);
numberOfServedUsersPerIteration = zeros(timeSteps,1);
minDistanceBetweenDronesPerIteration = zeros(timeSteps,1);

numberOfUsers = size(x,1);
%  m_initial = ceil(size(x,1)/maxUsersPerCluster)+1;
root= mcts_node_rewards_time(x,sigma,softmaxScale,maxUsersPerCluster,200,step,sinr_threshold,dronePower,numberOfUsers, desiredLinkReliability,0);

%%%%%%%%%%%%%%%
randomMu = unifrnd(s_input.V_POSITION_X_INTERVAL(1),s_input.V_POSITION_X_INTERVAL(2),...
    parameters.numberOfDrones,2);
randomSigma = 1.*dronePower.*ones(1,m_initial);
energyEfficiencyPerIterationRandom = zeros(timeSteps,1);
numberOfServedUsersPerIterationRandom = zeros(timeSteps,1);
time =1;
for timeIndex = timeIndexInitital+1:1:length(v_t)-1
    time = time +1;
    time_now = time*time_step
    if ~isempty(s_mobility_group)
        x=[reshape(leadersM(:,timeIndex,:),2,s_input.NB_NODES) reshape(membersM(:,timeIndex,:),2,totalMembers)]';
    else
        x=[reshape(leadersM(:,timeIndex,:),2,s_input.NB_NODES)]';
    end
    %     scatter(x(:,1),x(:,2))
    %     x(all(x==0,2),:)=[];
    if time_now == 1 || mod(time_now, parameters.mctsUpdateInterval) ==0
        maxScore_K_05_N_5 =0;
        reliabilities_K_05_N_5 = [];
        maxAllocationScore_K_05_N_5 = 0;
        maxEfficiencyScore_K_05_N_5 = 0;
        idxesToDelete = zeros(1,numel(childHistory_K_05_N_5));
        for idx =1:numel(childHistory_K_05_N_5)
            root.restart = 1;
            childHistory_K_05_N_5(idx).restart = 1;
            childHistory_K_05_N_5(idx).visits = 1;
            if childHistory_K_05_N_5(idx).visits < parameters.childVisistsThreshold
                idxesToDelete(idx) = idx;
            end
        end
        idxesToDelete(idxesToDelete==0)=[];
        for idxToDel = flip(idxesToDelete)
            childHistory_K_05_N_5(idxToDel) = [];
            
            statesHistory_K_05_N_5(idxToDel+1,:) = [];
            
        end
%         root= mcts_node_rewards_time(x,sigma,softmaxScale,maxUsersPerCluster,200,step,sinr_threshold,dronePower,numberOfUsers, desiredLinkReliability,0);
%         childHistory_K_05_N_5=[];
%         statesHistory_K_05_N_5 = zeros(1,m_initial);
        firstVisits = firstVisitsCounter;
        for iter=1:parameters.maxIterMcts
            root.update(x,root.merit);
%             while firstVisits == firstVisitsCounter
%                 root.update(x,root.merit);
%             end
%             firstVisits = firstVisitsCounter;
        end
    end
    
    maxSigma(maxSigma==0)=[];
    SINR = getSINRFinal(x,maxMu,maxSigma./1000,70.*ones(1,length(maxSigma))); %%%%Height is also inside EM_PULL
    [maxSINR, indices] = max(SINR,[],2);
    usersExceeding = sum(((maxSINR > 10.^(sinr_threshold./10)).*indices)==(1:m_initial))...
        - maxUsersPerCluster;
    usersExceeding(usersExceeding<0)=0;
    numberOfServedUsers = (sum(maxSINR > 10.^(sinr_threshold./10))- sum(usersExceeding));
    
    
    numberOfStatesVisitedMcts(timeIndex) = length(statesHistory_K_05_N_5);
    reliabilitiesPerIterationMcts(:,timeIndex) =sum(10.*log10(maxSINR)>sinrThresholds)/n;
    energyEfficiencyPerIteration(timeIndex) = (bandwidthPerUser.*sum(log2(maxSINR + 1)))./sum(maxSigma);
    spectralEfficiencyPerIteration(timeIndex) = mean(log2(maxSINR + 1));
    numberOfServedUsersPerIteration(timeIndex) = numberOfServedUsers;
    averageDronePower(timeIndex) = sum(maxSigma)/m_initial;
    minDistanceBetweenDronesPerIteration(timeIndex) = getMinDistance(maxMu);
    
    
    SinrRandom = getSINRFinal(x,randomMu,randomSigma./1000,70.*ones(1,length(randomSigma))); %%%%Height is also inside EM_PULL
    [maxSinrRandom, indicesRandom] = max(SinrRandom,[],2);
    usersExceedingRandom = sum(((maxSinrRandom > 10.^(sinr_threshold./10)).*indicesRandom)==(1:m_initial))...
        - maxUsersPerCluster;
    usersExceedingRandom(usersExceedingRandom<0)=0;
    numberOfServedUsersRandom = (sum(maxSinrRandom > 10.^(sinr_threshold./10))- sum(usersExceedingRandom));
    energyEfficiencyPerIterationRandom(timeIndex) = (bandwidthPerUser.*sum(log2(maxSinrRandom + 1)))./sum(randomSigma);
    numberOfServedUsersPerIterationRandom(timeIndex) = numberOfServedUsersRandom;
end