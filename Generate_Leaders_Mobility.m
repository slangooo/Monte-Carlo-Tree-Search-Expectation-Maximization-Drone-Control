function s_mobility = Generate_Mobility_MY(s_input)
global s_mobility_tmp;
global nodeIndex_tmp;
s_mobility.NB_NODES = s_input.NB_NODES;
s_mobility.SIMULATION_TIME = s_input.SIMULATION_TIME;

%TODO PREALLOC !!
for nodeIndex_tmp = 1:s_mobility.NB_NODES

    s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME = [];
    s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_X = [];
    s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_Y = [];
    
    X = unifrnd(s_input.V_POSITION_X_INTERVAL(1),s_input.V_POSITION_X_INTERVAL(2));
    Y = unifrnd(s_input.V_POSITION_Y_INTERVAL(1),s_input.V_POSITION_Y_INTERVAL(2));
    % Duration = 0;
    Time = 0;
    
    s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end+1) = 0;
    s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_X(end+1) = X;
    s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_Y(end+1) = Y;
    
    Walk_Leader(X,Y,Time,s_input);
    
    while (s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end) < s_input.SIMULATION_TIME-100*eps)
        Walk_Leader(s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_X(end),s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_Y(end),s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end),s_input);
    end  
end
s_mobility=s_mobility_tmp;
end

function Walk_Leader(x_tmp,y_tmp,Time,s_input)
global s_mobility_tmp;
global nodeIndex_tmp;
time_tmp = Time;
speed = unifrnd(s_input.V_SPEED_INTERVAL(1),s_input.V_SPEED_INTERVAL(2));
epsilon = 0.1;
if speed <epsilon
    duration_tmp = Out_adjustDuration_random_waypoint(time_tmp,unifrnd(s_input.V_PAUSE_INTERVAL(1),s_input.V_PAUSE_INTERVAL(2)),s_input);
    s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end+1) = time_tmp + duration_tmp;
    s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_X(end+1) =  x_tmp;
    s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_Y(end+1) =  y_tmp;
else
    duration_tmp = Out_adjustDuration_random_waypoint(time_tmp,unifrnd(s_input.V_WALK_INTERVAL(1),s_input.V_WALK_INTERVAL(2)),s_input);
    direction_tmp = unifrnd(s_input.V_DIRECTION_INTERVAL(1),s_input.V_DIRECTION_INTERVAL(2));
    distance_tmp = speed*duration_tmp;
    flag_mobility_finished = false;
    while (~flag_mobility_finished)
        x_dest = x_tmp + distance_tmp*cosd(direction_tmp);
        y_dest = y_tmp + distance_tmp*sind(direction_tmp);
        flag_mobility_was_outside = false;
        if (x_dest > s_input.V_POSITION_X_INTERVAL(2))
            flag_mobility_was_outside = true;
            new_direction = 180 - direction_tmp;
            x_dest = s_input.V_POSITION_X_INTERVAL(2);
            y_dest = y_tmp + diff([x_tmp x_dest])*tand(direction_tmp);
        end
        if (x_dest < s_input.V_POSITION_X_INTERVAL(1))
            flag_mobility_was_outside = true;
            new_direction = 180 - direction_tmp;
            x_dest = s_input.V_POSITION_X_INTERVAL(1);
            y_dest = y_tmp + diff([x_tmp x_dest])*tand(direction_tmp);
        end
        if (y_dest > s_input.V_POSITION_Y_INTERVAL(2))
            flag_mobility_was_outside = true;
            new_direction = -direction_tmp;
            y_dest = s_input.V_POSITION_Y_INTERVAL(2);
            x_dest = x_tmp + diff([y_tmp y_dest])/tand(direction_tmp);
        end
        if (y_dest < s_input.V_POSITION_Y_INTERVAL(1))
            flag_mobility_was_outside = true;
            new_direction = -direction_tmp;
            y_dest = s_input.V_POSITION_Y_INTERVAL(1);
            x_dest = x_tmp + diff([y_tmp y_dest])/tand(direction_tmp);
        end
        current_distance = abs(diff([x_tmp x_dest]) + 1i*diff([y_tmp y_dest]));
        current_duration = Out_adjustDuration_random_waypoint(time_tmp,current_distance/speed,s_input);
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_TIME(end+1) = time_tmp+current_duration;
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_X(end+1) = x_dest;
        s_mobility_tmp.VS_NODE(nodeIndex_tmp).V_POSITION_Y(end+1) = y_dest;
        if(flag_mobility_was_outside)
            time_tmp = time_tmp + current_duration;
            duration_tmp = duration_tmp - current_duration;
            distance_tmp = distance_tmp - current_distance;
            x_tmp = x_dest;
            y_tmp = y_dest;
            direction_tmp = new_direction;
        else
            flag_mobility_finished = true;
        end
    end
end
end

function duration = Out_adjustDuration_random_waypoint(time,duration,s_input)
if ((time+duration) >= s_input.SIMULATION_TIME)
    duration = s_input.SIMULATION_TIME - time;
end
end