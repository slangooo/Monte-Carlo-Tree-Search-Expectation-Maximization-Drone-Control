function mobility_group = Generate_Members_Mobility_MY(s_input,s_mobility)
global s_mobility_group;
global nodeIndex_leader;
global nodeIndex_group;

s_input.SIMULATION_TIME = s_input.SIMULATION_TIME;
for nodeIndex_leader = 1:s_input.NB_NODES
    s_mobility_group(nodeIndex_leader).NB_NODES_group = randi(s_input.MEMBERS_NB_INTERVAL);
    %s_mobility_group(nodeIndex_leader).RMAX= unifrnd(s_input.MEMBERS_RMAX_INTERVAL(1),s_input.MEMBERS_RMAX_INTERVAL(2));
    for nodeIndex_group=1:s_mobility_group(nodeIndex_leader).NB_NODES_group
        %INITIALIZE
        
        s_mobility_group(nodeIndex_leader).VS_NODE(nodeIndex_group).V_TIME = 0;
        %theta=unifrnd(s_input.MEMBERS_DIRECTION_INTERVAL(1), s_input.MEMBERS_DIRECTION_INTERVAL(2));
        %r= unifrnd(0,s_mobility_group(nodeIndex_leader).RMAX);
        time_in=0;
        %leader_x = interp1(s_mobility.VS_NODE(nodeIndex_leader).V_TIME, s_mobility.VS_NODE(nodeIndex_leader).V_POSITION_X,time_tmp);
        %leader_y = interp1(s_mobility.VS_NODE(nodeIndex_leader).V_TIME, s_mobility.VS_NODE(nodeIndex_leader).V_POSITION_Y,time_tmp);
        %s_mobility_group(nodeIndex_leader).VS_NODE(nodeIndex_group).V_POSITION_X = boundV(r*cos(theta)+leader_x,s_input.POSITION_X_INTERVAL);
        %s_mobility_group(nodeIndex_leader).VS_NODE(nodeIndex_group).V_POSITION_Y = boundV(r*sin(theta)+leader_y,s_input.POSITION_Y_INTERVAL);
        leader_x_in=s_mobility.VS_NODE(nodeIndex_leader).V_POSITION_X(1);
        leader_y_in=s_mobility.VS_NODE(nodeIndex_leader).V_POSITION_Y(1);
        s_mobility_group(nodeIndex_leader).VS_NODE(nodeIndex_group).V_POSITION_X = leader_x_in;
        s_mobility_group(nodeIndex_leader).VS_NODE(nodeIndex_group).V_POSITION_Y = leader_y_in;
        
        %         member_x_in=s_mobility_group(nodeIndex_leader).VS_NODE(nodeIndex_group).V_POSITION_X(end);
        %         member_y_in=s_mobility_group(nodeIndex_leader).VS_NODE(nodeIndex_group).V_POSITION_Y(end);
        [leader_x_in,leader_y_in]=Walk_Member(leader_x_in,leader_y_in,time_in,s_input,s_mobility);
        
        while (time_in < s_input.SIMULATION_TIME - eps*1e2)
            [leader_x_in,leader_y_in]= Walk_Member(leader_x_in,leader_y_in,time_in,s_input,s_mobility);
            time_in=s_mobility_group(nodeIndex_leader).VS_NODE(nodeIndex_group).V_TIME(end);
        end
        
    end
    
end
mobility_group=s_mobility_group;
end
function [leader_x_in,leader_y_in]=Walk_Member(leader_x_in,leader_y_in,time_in,s_input,s_mobility)
global s_mobility_group;
global nodeIndex_leader;
global nodeIndex_group;
%duration = Out_adjustDuration_random_waypoint(time_in,unifrnd(s_input.MEMBERS_WALK_INTERVAL(1),s_input.MEMBERS_WALK_INTERVAL(2)),s_input);
%time_final=time_in+duration;

if numel(s_mobility_group(nodeIndex_leader).VS_NODE(nodeIndex_group).V_TIME)+1 > numel(s_mobility.VS_NODE(nodeIndex_leader).V_TIME)
    time_final=s_mobility.VS_NODE(nodeIndex_leader).V_TIME(end);
else
    time_final=s_mobility.VS_NODE(nodeIndex_leader).V_TIME(numel(s_mobility_group(nodeIndex_leader).VS_NODE(nodeIndex_group).V_TIME)+1);
end

while s_mobility.VS_NODE(nodeIndex_leader).V_TIME(end) - s_mobility.VS_NODE(nodeIndex_leader).V_TIME(end-1)<0.1
    s_mobility.VS_NODE(nodeIndex_leader).V_TIME(end)=[];
    s_mobility.VS_NODE(nodeIndex_leader).V_POSITION_X(end)=[];
    s_mobility.VS_NODE(nodeIndex_leader).V_POSITION_Y(end)=[];
end
try
    leader_x_final= interp1(s_mobility.VS_NODE(nodeIndex_leader).V_TIME, s_mobility.VS_NODE(nodeIndex_leader).V_POSITION_X,time_final);
    leader_y_final= interp1(s_mobility.VS_NODE(nodeIndex_leader).V_TIME, s_mobility.VS_NODE(nodeIndex_leader).V_POSITION_Y,time_final);
catch
    bla = 5;
end
speed_deviation= norm([leader_x_final leader_y_final]-[leader_x_in,leader_y_in])+unifrnd(-1,1)*s_input.SDR*s_input.V_SPEED_INTERVAL(2);
angle_deviation= getAngle([leader_x_final leader_y_final],[leader_x_in,leader_y_in])+unifrnd(-1,1)*s_input.ADR*s_input.V_DIRECTION_INTERVAL(2);

s_mobility_group(nodeIndex_leader).VS_NODE(nodeIndex_group).V_POSITION_X(end+1)= boundV(leader_x_in + speed_deviation*cosd(angle_deviation),s_input.V_POSITION_X_INTERVAL);
s_mobility_group(nodeIndex_leader).VS_NODE(nodeIndex_group).V_POSITION_Y(end+1)= boundV(leader_y_in + speed_deviation*sind(angle_deviation),s_input.V_POSITION_Y_INTERVAL);
s_mobility_group(nodeIndex_leader).VS_NODE(nodeIndex_group).V_TIME(end+1)=time_final;

leader_x_in=leader_x_final;
leader_y_in=leader_y_final;

end

function duration = Out_adjustDuration_random_waypoint(time,duration,s_input)
if ((time+duration) >= s_input.SIMULATION_TIME)
    duration = s_input.SIMULATION_TIME - time;
end
end