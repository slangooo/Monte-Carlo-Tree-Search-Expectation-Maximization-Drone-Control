classdef mcts_node_rewards_K_05 < handle
    % A node for a tree in a MCTS
    %K = 0.01
    properties
        K=0.01;
        discountFactor =1;
        maxSigma = 500;
        reliability_score = 0;
        allocation_score=0;
        efficiency_score =0;
        merit = 0;
        maxIter = 0;
        softmaxScale = 0;
        sigma=[];
        maxUsersPerCluster = [];
        m=0;
        step=0;
        childs = {};
        num_childs = 0;
        x=[];
        visits =0;
        scores =[];
        restart =0;
        sinr_threshold=0;
        maxSINR = [];
        numberOfServedUsers = 0;
        max_allocation_score = 0;
        max_efficiency_score = 0;
        EmIter = 0;
        numberOfUsers = 0;
        desiredLinkReliability = 1;
        
    end
    methods
        function obj = mcts_node_rewards_K_05(x,sigma,softmaxScale,maxUsersPerCluster,maxIter,step,sinr_threshold,maxSigma,numberOfUsers,desiredLinkReliability)
            % for a new node simulate a score using MC
            obj.desiredLinkReliability = desiredLinkReliability;
            obj.numberOfUsers = numberOfUsers;
            obj.maxSigma=maxSigma;
            obj.sinr_threshold = sinr_threshold;
            obj.step = step;
            obj.m = length(sigma(sigma>0));
            obj.sigma = sigma;
            local_sigma =sigma(sigma>0);
            obj.maxIter = maxIter;
            obj.max_allocation_score = mcts_node_rewards_K_05.getMaxAllocationScore();
            obj.max_efficiency_score = mcts_node_rewards_K_05.getMaxEfficiencyScore();
            obj.maxUsersPerCluster = maxUsersPerCluster;
            obj.scores = zeros(1,obj.m);
            
            %Invalid state reached
            if any(obj.sigma<0) || any(obj.sigma > obj.maxSigma)
                obj.merit = -0.5;
                return
            end
            
            %Perform EM
            [obj.efficiency_score,obj.allocation_score,obj.reliability_score, obj.maxSINR, obj.numberOfServedUsers,obj.EmIter]...
                = EM_PULL_CNST_SIGMA1 (x,obj.m,local_sigma,softmaxScale,maxUsersPerCluster,maxIter,sinr_threshold);
            
            %Update the average number of EM iterations
            obj.UpdateEmIter();
            
            %Obtain node figure of merit
            obj.merit= obj.reliability_score/ceil(sqrt(obj.allocation_score+1));
            
            %If link reliability is achieved improve energy efficiency
            if obj.numberOfServedUsers/obj.numberOfUsers > obj.desiredLinkReliability
                obj.merit = obj.merit + obj.efficiency_score;
            end
            
            % %             %Check for win {more users served / better efficiency AND same
            % %             %number of users served AND reliability achieved}
            % %             if obj.numberOfServedUsers>obj.max_allocation_score ||...
            % %                     (obj.efficiency_score > obj.max_efficiency_score ...
            % %                     && obj.numberOfServedUsers == obj.max_allocation_score...
            % %                     && obj.reliability_score - 1 > obj.desiredLinkReliability)
            % %                 %Update maximum number of users served
            % %                 obj.setMaxAllocationScore();
            % %                 %Update transmission powers that achieved max state
            % %                 obj.setMaxState();
            % %                 %Update SINRs of max state and maximum reliability score
            % %                 obj.setMaxScore();
            % %             end
            %Cell to hold future childs
            obj.childs = cell(obj.m,1);
        end
        
        function [child_merit,iter] = update(obj,x)
            obj.max_allocation_score = mcts_node_rewards_K_05.getMaxAllocationScore();
            obj.max_efficiency_score = mcts_node_rewards_K_05.getMaxEfficiencyScore();
            if obj.restart && ~(any(obj.sigma<0) || any(obj.sigma > obj.maxSigma))
                obj.restart =0;
                [obj.efficiency_score,obj.allocation_score,obj.reliability_score, obj.maxSINR, obj.numberOfServedUsers,obj.EmIter] = EM_PULL_CNST_SIGMA1 (x,obj.m,obj.local_sigma,obj.softmaxScale,obj.maxUsersPerCluster,obj.maxIter);
                %Update the average number of EM iterations
                obj.UpdateEmIter();
                
                %Obtain node figure of merit
                obj.merit= obj.reliability_score/ceil(sqrt(obj.allocation_score+1));
                
                %If link reliability is achieved improve energy efficiency
                if obj.numberOfServedUsers/obj.numberOfUsers > obj.desiredLinkReliability
                    obj.merit = obj.merit + obj.efficiency_score;
                end
                
                
% %                 %Check for win {more users served / better efficiency AND same
% %                 %number of users served AND reliability achieved}
% %                 if obj.numberOfServedUsers>obj.max_allocation_score ||...
% %                         (obj.efficiency_score > obj.max_efficiency_score ...
% %                         && obj.numberOfServedUsers == obj.max_allocation_score...
% %                         && obj.reliability_score - 1 > obj.desiredLinkReliability)
% %                     %Update maximum number of users served
% %                     obj.setMaxAllocationScore();
% %                     %Update transmission powers that achieved max state
% %                     obj.setMaxState();
% %                     %Update SINRs of max state and maximum reliability score
% %                     obj.setMaxScore();
% %                 end
            end
            
            %Increment visits counter
            obj.visits = 1+ obj.visits;
            
            %Check for loss
            if (obj.visits ==1 || obj.restart) &&...
                    obj.numberOfServedUsers<obj.max_allocation_score ||...
                    any(obj.sigma<0)|| any(obj.sigma > obj.maxSigma)
                child_merit = 0;
                obj.restart =0;
                iter =0;
                return
            end
            
            %Check for win
            if (obj.visits ==1 || obj.restart) &&...
                    (obj.numberOfServedUsers>obj.max_allocation_score ||...
                    (obj.efficiency_score > obj.max_efficiency_score ...
                    && obj.numberOfServedUsers == obj.max_allocation_score...
                    && obj.reliability_score - 1 > obj.desiredLinkReliability))
                iter =0;
                child_merit = obj.merit;
                obj.restart =0;
                %Update maximum number of users served
                obj.setMaxAllocationScore();
                %Update transmission powers that achieved max state
                obj.setMaxState();
                %Update SINRs of max state and maximum reliability score
                obj.setMaxScore();
                return
            end
            
            % % % %             %If first visit
            % % % %             if (obj.visits ==1 || obj.restart) &&...
            % % % %                     (obj.allocation_score ==0 ||...
            % % % %                     obj.numberOfServedUsers>obj.max_allocation_score)
            % % % %                 iter =0;
            % % % %                 child_merit = obj.merit;
            % % % %                 if obj.numberOfServedUsers<obj.max_allocation_score
            % % % %                     child_merit = 0;
            % % % %                 end
            % % % %                 obj.restart =0;
            % % % %                 return
            % % % %             elseif
            % % % %                 obj.restart =0;
            % % % %                 iter =0;
            % % % %                 child_merit = 0;
            % % % %                 return
            % % % %             end
            obj.restart =0;
            
            %Make sure all children nodes are created
            if obj.num_childs == 0
                for i = 1:size(obj.sigma,2)
                    obj.num_childs = obj.num_childs + 1;
                    sigma_to_explore = obj.sigma;
                    sigma_to_explore(i) = sigma_to_explore(i) + obj.step;
                    index =  obj.checkIfStateOccurred(sigma_to_explore);
                    if index ~= 0
                        obj.childs{obj.num_childs} = obj.getChildWithIndex(index);
                    else
                        obj.childs{obj.num_childs} = mcts_node_rewards_K_05(x,sigma_to_explore,obj.softmaxScale,...
                            obj.maxUsersPerCluster,obj.maxIter,obj.step,obj.sinr_threshold,obj.maxSigma, obj.numberOfUsers, obj.desiredLinkReliability);
                        obj.addChildToHistory(obj.childs{obj.num_childs});
                    end
                end
            end
            
            
            
            % explore a potential state 
            sigma_to_explore = obj.explore();
            logicalIndex = logical(sigma_to_explore - obj.sigma);
            
            try
                [child_merit,iter] = obj.childs{min(find(logicalIndex),obj.num_childs)}.update(x);
            catch
                child_merit = -0.5;
                iter = 0;
            end
            iter = iter +1;
            obj.update_score(child_merit,iter);
        end
        function Path = explore(obj)
            %select a next state to explore, that may or may not be visited
            for i = 1:obj.num_childs
                obj.scores(i)=obj.childs{i}.merit+...
                    obj.K*sqrt(log(obj.visits)/(obj.childs{i}.visits+1));
            end
            [~,index]=max(obj.scores);
            new_sigma = obj.sigma;
            new_sigma(index) = new_sigma(index) + obj.step;
            Path = new_sigma;
        end
        
        function bool = state_equals(obj, test_state)
            % returns true if the nodes state is equal to test_state
            bool=all(obj.sigma==test_state);
        end
        
        function update_score(obj, child_merit, iter)
            % obj.discountFactor = 0.9; %EDIT: discount factor
            if  child_merit > 0
                obj.merit = obj.merit + obj.discountFactor^iter*(child_merit-obj.merit)/obj.visits;
            end
        end
        
        function setMaxScore(obj)
            global maxScore_K_05_N_5;
            global reliabilities_K_05_N_5;
            maxScore_K_05_N_5=obj.reliability_score;
            reliabilities_K_05_N_5 = obj.maxSINR;
        end
        function UpdateEmIter(obj)
            global EmIter
            global EmIterCounter
            EmIterCounter = EmIterCounter + 1;
            EmIter = EmIter + (obj.EmIter - EmIter)/EmIterCounter;
        end
        function setMaxState(obj)
            global maxSigma;
            global maxM;
            maxSigma=obj.sigma;
            maxM=obj.m;
        end
        function setMaxAllocationScore(obj)
            global maxAllocationScore_K_05_N_5;
            global maxEfficiencyScore_K_05_N_5;
            maxAllocationScore_K_05_N_5=obj.numberOfServedUsers;
            maxEfficiencyScore_K_05_N_5=obj.efficiency_score;
        end
    end
    methods(Static)
        function score = getMaxAllocationScore()
            global maxAllocationScore_K_05_N_5;
            score = maxAllocationScore_K_05_N_5;
        end
        function score = getMaxEfficiencyScore()
            global maxEfficiencyScore_K_05_N_5;
            score = maxEfficiencyScore_K_05_N_5;
        end
        function score = getMaxScore()
            global maxScore_K_05_N_5;
            score = maxScore_K_05_N_5;
        end
        function index =checkIfStateOccurred(state)
            global statesHistory_K_05_N_5;
            [bool, index] = ismember(state,statesHistory_K_05_N_5,'rows');
            if ~bool
                statesHistory_K_05_N_5 = [statesHistory_K_05_N_5; state];
            end
            index = max(0,index-1);
        end
        function child = getChildWithIndex(index)
            global childHistory_K_05_N_5;
            child = childHistory_K_05_N_5(index);
        end
        function addChildToHistory(child)
            global childHistory_K_05_N_5;
            childHistory_K_05_N_5 = [childHistory_K_05_N_5 child];
        end
    end
end