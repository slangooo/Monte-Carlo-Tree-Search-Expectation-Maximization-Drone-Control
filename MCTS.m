classdef mcts_node < handle
    % A node for a tree in a MCTS
    properties
        state = []; %some state of the search space that identifies the node
        reliability_score = 0;
        merits_updated_from_child = [];
        allocation_score=[];
        merit_scores = [];
        merit = 0;
        maxIter = 0;
        softmaxScale = 0;
        sigma=[];
        maxUsersPerCluster = [];
        m=0;
        childs = {};
        num_childs = 0;
        max_reliability_score=[];
        x=[];
    end
    methods
        function obj = node(x,sigma,softmaxScale,maxUsersPerCluster,maxIter, parent_merits)
            % for a new node simulate a score using MC
            obj.m = length(sigma(sigma~=0));
            obj.merits_updated_from_child = zeros(obj.m,1);
            obj.x=x;
            obj.sigma = sigma;
            obj.maxIter = maxIter;
            obj.softmaxScale = softmaxScale;
            obj.merit_scores = parent_merits;
            obj.max_reliability_score=node.getMaxScore();
            obj.maxUsersPerCluster = maxUsersPerCluster;
            %                                                                             (x,m,sigma,softmaxScale,maxUsersPerCluster,maxIter)
            [obj.allocation_score,obj.reliability_score] = EM_PULL_CNST_SIGMA (x,obj.m,sigma,softmaxScale,maxUsersPerCluster,maxIter);
            %             [obj.mu,obj.sigma,obj.m,obj.Pj,~,obj.I_score,obj.off,obj.score]= EM_Monte4(x,sigma,mu,Pj,m,scaleSigma); % TODO scaleSigma from State;
            obj.merit= obj.reliability_score/ceil(log(obj.allocation_score+exp(1)));
            if  obj.reliability_score>obj.max_reliability_score && obj.allocation_score == 0
                obj.setMaxScore();
                obj.setMaxState();
            end
            obj.childs = cell(obj.m,1);
        end
        
        function value= update(obj)
            % explore a potential state (add new drone)
            sigma_to_explore = obj.explore();
            
            %check if state has already been visited
            terminate = false;
            idx = 1;
            while idx <= obj.num_childs && ~terminate
                if obj.childs{idx}.state_equals(sigma_to_explore)
                    terminate = true;
                end
                idx = idx + 1;
            end
            %preform the according action based on search
            if ~terminate
                % state has never been visited
                % this action terminates the update recursion
                % and creates a new leaf
                obj.num_childs = obj.num_childs + 1;
                obj.childs{obj.num_childs} = node(obj.x,sigma_to_explore,obj.softmaxScale,...
                    obj.maxUsersPerCluster,obj.maxIter,obj.merit_scores);
                child_merit = obj.childs{obj.num_childs}.merit;
%                 child_allocation_score =  obj.childs{obj.num_childs}.allocation_score;
                obj.update_score(child_merit,sigma_to_explore);
            else
                % state has been visited at least once
                value = obj.childs{idx-1}.update();
                obj.update_score(value);
            end
        end
        function Path = explore(obj)
            %select a next state to explore, that may or may not be visited
            Path=pickPath(softmax(obj.merit_scores/softmax_scale));
        end
        
        function bool = do_exploration(obj)
            % decide if this node should be explored or exploited
            tmp = obj.allocation_score / obj.maxUsersPerCluster;
            if tmp>=1
                bool = 1;
            else
                bool = pickPath([tmp 1]) - 1;
            end
        end
        
        function bool = state_equals(obj, test_state)
            % returns true if the nodes state is equal to test_state
            bool=~any(~(obj.state==test_state));
        end
        function update_score(obj, child_merit, sigma_explored)
            index = find(logical(obj.sigma - sigma_explored));
            if obj.merit > child_merit
                obj.merit = child_merit;
            end
            if ~obj.merits_updated_from_child(index)
                obj.merit_scores(index) = child_merit;
                obj.merits_updated_from_child(index)=1;
            else
                if obj.merit_scores(index)> child_merit
                    obj.merit_scores(index) = child_merit;
                end
            end            
        end
        function setMaxScore(obj)
            global maxScore;
            maxScore=obj.reliability_score;
        end
        function setMaxState(obj)
            global maxSigma;
            global maxM;
            maxSigma=obj.sigma;
            maxM=obj.m;
        end
    end
    methods(Static)
        function score = getMaxScore()
            global maxScore;
            score = maxScore;
        end
    end
end