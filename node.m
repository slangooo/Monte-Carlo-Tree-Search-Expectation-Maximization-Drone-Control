classdef node < handle
    % A node for a tree in a MCTS
    properties
        state = []; %some state of the search space that identifies the node
        reliability_score = 0;
        m_max = 0;
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
        function obj = node(x,m_max,sigma,softmaxScale,maxUsersPerCluster,maxIter, parent_merits)
            % for a new node simulate a score using MC
            obj.x=x;
            obj.sigma = sigma;
            obj.m_max = m_max;
            obj.maxIter = maxIter;
            obj.softmaxScale = softmaxScale;
            obj.merit_scores = parent_merits;
            obj.max_reliability_score=node.getMaxScore();
            obj.maxUsersPerCluster = maxUsersPerCluster;
%                                                                             (x,m,sigma,softmaxScale,maxUsersPerCluster,maxIter)
            obj.m = length(sigma(sigma~=0));
            [obj.allocation_score,obj.reliability_score] = EM_PULL_CNST_SIGMA (x,obj.m,sigma,softmaxScale,maxUsersPerCluster,maxIter);
%             [obj.mu,obj.sigma,obj.m,obj.Pj,~,obj.I_score,obj.off,obj.score]= EM_Monte4(x,sigma,mu,Pj,m,scaleSigma); % TODO scaleSigma from State;

            if  obj.reliability_score>obj.max_reliability_score && obj.allocation_score == 0
                obj.setMaxScore();
                obj.setMaxState();
            end
            obj.childs = cell(m_max,1);
        end
        function value = update(obj)
            
            % update the this node using MC recursively
            if obj.do_exploration() || obj.num_childs == 0
                % explore a potential state (add new drone)
                state_to_explore = obj.explore();
                
                %check if state has already been visited
                terminate = false;
                idx = 1;
                while idx <= obj.num_childs && ~terminate
                    if obj.childs{idx}.state_equals(state_to_explore)
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
                    obj.childs{obj.num_childs} = node(obj.x,obj.m_max,state_to_explore,obj.softmaxScale,...
                        obj.maxUsersPerCluster,obj.maxIter,obj.merit_scores);
                    child_merit = obj.childs{obj.num_childs}.merit;
                    obj.update_score(child_merit);
                else
                    % state has been visited at least once
                    value = obj.childs{idx-1}.update();
                    obj.update_score(value);
                end
            else
                % exploit what we know already
                dist=zeros(obj.num_childs);
                for idx = 1:obj.num_childs
                    dist(idx)=obj.childs{idx}.score;
                end
                choice=softmax(((dist-max(dist))/(1*min(dist))));%%%%%%%%%
                choice=pickPath(choice);
                value = obj.childs{choice}.update();
                obj.update_score(value);
            end
            value = obj.score;
        end
        function Path = explore(obj)
            %select a next state to explore, that may or may not be visited
            Path=pickPath(obj.merit_scores);
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
        function update_score(obj, value)
            % updates the score based on some value
            if value>obj.score
                obj.score=value;
            end
        end
        function TurnOff (obj)
            jjj=0;
            for j=1:obj.m
                jj=j-jjj;
                if(obj.sigma(1,1,jj)<obj.MinSigma)
                    obj.m=obj.m-1;
                    obj.mu(jj,:)=[];
                    obj.sigma(:,:,jj)=[];
                    obj.Pj(jj)=[];
                    jjj=jjj+1;
                    obj.state(obj.state==jj)=0;
                    obj.state(obj.state>jj)=obj.state(obj.state>jj)-1;
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