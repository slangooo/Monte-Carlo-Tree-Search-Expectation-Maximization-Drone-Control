% clc
% clear
% maxUsersPerCluster_local=[250 350 450];
% SINR_Thresh_local = [4 2 0];
% % maxUsersPerCluster_local=[450 250];
% % SINR_Thresh_local = [4];
% h=10;
% pathLossExponent=2;
% refPower = -30; %reference power at 1m dbm
% Noise= -100; % receiver noise power dBm;
% sigma_max=1500; % mwatts 
% sigma_step=-150; %mwatts 
% m_initial = 5;
% ii=0;
% sinr_achieved=zeros(500,10,4,numel(SINR_Thresh_local),numel(SINR_Thresh_local)); %SECONDS - ITERATION - MCTS_MODES
% number_of_states_visited = zeros(500,10,4,numel(SINR_Thresh_local),numel(SINR_Thresh_local));%SECONDS - ITERATION - MCTS_MODES

while ii<1000
    ii= ii+1
    vars = {'s_mobility_group', 's_mobility'};
    clear(vars{:});
    s_input = struct('V_POSITION_X_INTERVAL',[0 500],...%(m)
        'V_POSITION_Y_INTERVAL',[0 500],...%(m)
        'V_SPEED_INTERVAL',[0.8 1.4],...%(m/s)
        'V_PAUSE_INTERVAL',[0 1],...%pause time (s)
        'V_WALK_INTERVAL',[2.00 6.00],...%walk time (s)
        'V_DIRECTION_INTERVAL',[0 360],...%(degrees)
        'SIMULATION_TIME',60,...%(s)
        'NB_NODES',randi([8 10]),...
        'MEMBERS_NB_INTERVAL', [90 100],...
        'MEMBERS_WALK_INTERVAL', [1 3],...
        'SDR', 20,...  %speed member deviation ratio
        'ADR',1); %angle member devation ratio
    
    %
    s_mobility = Generate_Leaders_Mobility(s_input);
    s_mobility_group = Generate_Members_Mobility_Synchronous(s_input,s_mobility);
    time_step = 0.5;%(s)

    for sinr_index=1:numel(SINR_Thresh_local)
        sinr_threshold=SINR_Thresh_local(sinr_index)
        parfor (max_users_index=1:numel(maxUsersPerCluster_local))
             maxUsersPerCluster=maxUsersPerCluster_local(max_users_index)
             [number_of_states_visited(:,:,:,sinr_index,max_users_index),sinr_achieved(:,:,:,sinr_index,max_users_index)...
                 ]=test_Simulate...
                 (s_mobility,s_input,time_step,s_mobility_group,sigma_max,m_initial,...
                 sigma_step, maxUsersPerCluster, sinr_threshold)
        end
    end
    if ii==1
        number_of_states_visited_total = number_of_states_visited;
        sinr_achieved_total = sinr_achieved;
    else
        number_of_states_visited_total = number_of_states_visited_total +...
            (number_of_states_visited - number_of_states_visited_total)./ii;
        sinr_achieved_total = sinr_achieved_total + (sinr_achieved - sinr_achieved_total)./ii;
    end
    save('simulation_results','number_of_states_visited_total','sinr_achieved_total');
end