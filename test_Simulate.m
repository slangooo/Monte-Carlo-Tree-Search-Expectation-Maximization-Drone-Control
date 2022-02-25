function [number_of_states_visited,sinr_achieved]=test_Simulate(s_mobility,s_input,time_step,s_mobility_group,sigma_max,m_initial,sigma_step,maxUsersPerCluster,sinr_threshold)
% global maxUsersPerCluster
% global SINR_Thresh
v_t = 0:time_step:s_input.SIMULATION_TIME;
store= struct('v_x', zeros(1,1+s_input.SIMULATION_TIME/time_step),...
    'v_y',zeros(1,1+s_input.SIMULATION_TIME/time_step));
leaders=repmat(store,1,s_input.NB_NODES);
cellGroup=struct2cell(s_mobility_group);
Nmembers=cell2mat(cellGroup(1,:));
cellGroup=[];
members=repmat(store,1,sum(Nmembers));
countMembers=0;
for leaderIndex = 1:s_input.NB_NODES
    leaders(leaderIndex).v_x = interp1(s_mobility.VS_NODE(leaderIndex).V_TIME,s_mobility.VS_NODE(leaderIndex).V_POSITION_X,v_t);
    leaders(leaderIndex).v_y = interp1(s_mobility.VS_NODE(leaderIndex).V_TIME,s_mobility.VS_NODE(leaderIndex).V_POSITION_Y,v_t);
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
leadersM=cell2mat(struct2cell(leaders));
membersM=cell2mat(struct2cell(members));
vars = {'cellGroup', 'store', 'countMembers','members','leaders'};
clear(vars{:});
totalMembers=sum(Nmembers);

sinr_achieved=zeros(500,10,4); %SECONDS - ITERATION - MCTS_MODES
number_of_states_visited = zeros(500,10,4);%SECONDS - ITERATION - MCTS_MODES

m = m_initial;
sigma = 1.*sigma_max.*ones(1,m);
softmaxScale = 50;
maxUsersPerCluster = maxUsersPerCluster;
sinr_threshold = sinr_threshold;
maxIter = 150;
step = sigma_step;

% global maxScore_K_05_N_5
% global statesHistory_K_05_N_5
% global childHistory_K_05_N_5
% childHistory_K_05_N_5=[];
% statesHistory_K_05_N_5 = zeros(1,m_initial);
% maxScore_K_05_N_5 =0;

global maxScore_K_09_N_5
global statesHistory_K_09_N_5
global childHistory_K_09_N_5
childHistory_K_09_N_5=[];
statesHistory_K_09_N_5 = zeros(1,m_initial);
maxScore_K_09_N_5 =0;

global maxScore_K_09_N_2
global statesHistory_K_09_N_2
global childHistory_K_09_N_2
childHistory_K_09_N_2=[];
statesHistory_K_09_N_2 = zeros(1,m_initial);
maxScore_K_09_N_2 =0;

global maxScore_K_05_N_2
global statesHistory_K_05_N_2
global childHistory_K_05_N_2
childHistory_K_05_N_2=[];
statesHistory_K_05_N_2 = zeros(1,m_initial);
maxScore_K_05_N_2 =0;


timeIndexInitital = 15;
x=[reshape(leadersM(:,timeIndexInitital,:),2,s_input.NB_NODES) reshape(membersM(:,timeIndexInitital,:),2,totalMembers)]';
x(all(x==0,2),:)=[];
time =1;

% root_05_05= mcts_node_rewards_K_05(x,sigma,softmaxScale,maxUsersPerCluster,maxIter,step,sinr_threshold,sigma_max);
root_09_05= mcts_node_rewards_K_09(x,sigma,softmaxScale,maxUsersPerCluster,maxIter,step,sinr_threshold,sigma_max);
root_05_02= mcts_node_rewards_K_05_N_02(x,sigma,softmaxScale,maxUsersPerCluster,maxIter,step,sinr_threshold,sigma_max);
root_09_02= mcts_node_rewards_K_09_N_02(x,sigma,softmaxScale,maxUsersPerCluster,maxIter,step,sinr_threshold,sigma_max);

% roots ={root_05_05;root_09_05;root_05_02;root_09_02};
roots ={root_09_05;root_05_02;root_09_02};
% roots ={root_05_02};
num_updates = 60;
update_step =num_updates/10;


for i = 1:num_updates
    for root_num =1:numel(roots)
        roots{root_num}.update(x);
        %         state_histories ={statesHistory_K_05_N_5;statesHistory_K_09_N_5;statesHistory_K_05_N_2;statesHistory_K_09_N_2};
        %         max_scores = {maxScore_K_05_N_5;maxScore_K_09_N_5;maxScore_K_05_N_2;maxScore_K_09_N_2};
        state_histories ={statesHistory_K_09_N_5;statesHistory_K_05_N_2;statesHistory_K_09_N_2};
        max_scores = {maxScore_K_09_N_5;maxScore_K_05_N_2;maxScore_K_09_N_2};
        if mod(i,update_step)==0
            sinr_achieved(time,i/update_step,root_num)=max_scores{root_num}-1;
            number_of_states_visited(time,i/update_step,root_num) = length(state_histories{root_num});
        end
    end
end
% child_histories= {childHistory_K_05_N_5;childHistory_K_09_N_5 ;childHistory_K_05_N_2 ; childHistory_K_05_N_2};
child_histories= {childHistory_K_09_N_5 ;childHistory_K_05_N_2 ; childHistory_K_05_N_2};
for root_num =1:numel(roots)
    for idx =1:numel(child_histories{root_num})
        child_histories{root_num}(idx).restart = 0;
        if root_num == 2 || root_num ==3
            child_histories{root_num}(idx).visits = max(1,ceil(child_histories{root_num}(idx).visits/2));
        else
            child_histories{root_num}(idx).visits = max(1,ceil(child_histories{root_num}(idx).visits/5));
        end
    end
end

num_updates = 10;
update_step =num_updates/10;

for timeIndex = timeIndexInitital+1:1:length(v_t)-3
    time = time +1
    x=[reshape(leadersM(:,timeIndex,:),2,s_input.NB_NODES) reshape(membersM(:,timeIndex,:),2,totalMembers)]';
    x(all(x==0,2),:)=[];
    for i = 1:num_updates
        for root_num =1:numel(roots)
            roots{root_num}.update(x);
            %             state_histories ={statesHistory_K_05_N_5;statesHistory_K_09_N_5;statesHistory_K_05_N_2;statesHistory_K_09_N_2};
            %             max_scores = {maxScore_K_05_N_5;maxScore_K_09_N_5;maxScore_K_05_N_2;maxScore_K_09_N_2};
            state_histories ={statesHistory_K_09_N_5;statesHistory_K_05_N_2;statesHistory_K_09_N_2};
            max_scores = {maxScore_K_09_N_5;maxScore_K_05_N_2;maxScore_K_09_N_2};
            
            if mod(i,update_step)==0
                sinr_achieved(time,i/update_step,root_num)=max_scores{root_num}-1;
                number_of_states_visited(time,i/update_step,root_num) = length(state_histories{root_num});
            end
        end
    end
    %     child_histories= {childHistory_K_05_N_5;childHistory_K_09_N_5 ;childHistory_K_05_N_2 ; childHistory_K_05_N_2};
    child_histories= {childHistory_K_09_N_5 ;childHistory_K_05_N_2 ; childHistory_K_05_N_2};
    for root_num =1:numel(roots)
        for idx =1:numel(child_histories{root_num})
            child_histories{root_num}(idx).restart = 0;
            if root_num == 2 || root_num == 3
                child_histories{root_num}(idx).visits = max(1,ceil(child_histories{root_num}(idx).visits/2));
            else
                child_histories{root_num}(idx).visits = max(1,ceil(child_histories{root_num}(idx).visits/5));
            end
        end
    end
end
end