function test_Animate(s_mobility,s_input,time_step,s_mobility_group)
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
        members(groupIndex+countMembers).v_x=interp1(s_mobility_group(leaderIndex).VS_NODE(groupIndex).V_TIME,s_mobility_group(leaderIndex).VS_NODE(groupIndex).V_POSITION_X,v_t);
        members(groupIndex+countMembers).v_y=interp1(s_mobility_group(leaderIndex).VS_NODE(groupIndex).V_TIME,s_mobility_group(leaderIndex).VS_NODE(groupIndex).V_POSITION_Y,v_t);
    end
    countMembers=countMembers+s_mobility_group(leaderIndex).NB_NODES_group;
end

leadersM=cell2mat(struct2cell(leaders));
membersM=cell2mat(struct2cell(members));
vars = {'cellGroup', 'store', 'countMembers','members','leaders'};
clear(vars{:});

figure;
hold on;
totalMembers=sum(Nmembers);
% vh_node_pos=zeros(1,s_input.NB_NODES);
% vh_member_pos =zeros(1,totalMembers);
% for leaderIndex = 1:s_input.NB_NODES
%     vh_node_pos(leaderIndex) = plot(leadersM(1,1,leaderIndex),leadersM(2,1,leaderIndex),'*','color',[0.3 0.3 1]);
% end
% for memberIndex=1:totalMembers
%     vh_member_pos(memberIndex)=plot(membersM(1,1,memberIndex),membersM(2,1,memberIndex),'o','color',[0.3 1 0.3]);
% end
% 
% title(cat(2,'Simulation time (sec): ',num2str(s_input.SIMULATION_TIME)));
% xlabel('X (meters)');
% ylabel('Y (meters)');
% title('Radom Waypoint mobility');
% ht = text(10,30,cat(2,'Time (sec) = 0'));
% axis([[10 30] [10 30]]);
% hold on;

% initEMFlag=1;
% initPlotFlag=1;
% m=6;
for timeIndex = 1:length(v_t)
    t = v_t(timeIndex);
%     set(ht,'String',cat(2,'Time (sec) = ',num2str(t,4)));
%     for leaderIndex = 1:s_input.NB_NODES
%         set(vh_node_pos(leaderIndex),'XData',leadersM(1,timeIndex,leaderIndex),'YData',leadersM(2,timeIndex,leaderIndex));
%     end
%     for memberIndex=1:totalMembers
%         set(vh_member_pos(memberIndex),'XData',membersM(1,timeIndex,memberIndex),'YData',membersM(2,timeIndex,memberIndex));
%     end
%     drawnow;
    if timeIndex*time_step>5
        x=[reshape(leadersM(:,timeIndex,:),2,s_input.NB_NODES) reshape(membersM(:,timeIndex,:),2,totalMembers)];
        
%         if initEMFlag
%             sigma=repmat([2 2].*eye(2),1,1,m);
%             mu=[];
%             Pj=[];
%             initEMFlag=0;
%         end
%         [mu,sigma,Allocs,Pj]=EM_FLAG(x',sigma,15,2,mu,Pj,m,4,0.00000001,0.1,0,2);
        
        
%         for j=1:m
%             aa=sigma(1,1,j); % horizontal radius
%             bb=sigma(2,2,j); % vertical radius
%             x0=mu(j,1); % x0,y0 ellipse centre coordinates
%             y0=mu(j,2);
%             tt=-pi:0.01:pi;
%             xp=x0+aa*cos(tt);
%             yp=y0+bb*sin(tt);
%             if initPlotFlag
%                 pp(j)=plot(xp,yp);
%                 if j==m
%                     initPlotFlag=0;
%                 end
%             else
%                 set(pp(j),'XData',xp,'YData',yp);
%             end
%         end
    end
    
end
end