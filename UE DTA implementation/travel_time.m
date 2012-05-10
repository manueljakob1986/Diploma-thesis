function [act_travel_time,Arriving,Departing] = travel_time(x_veh,flow,paths,discretesplits_nodes,dt,source,num_demands)

%Here we are going to compute the actual travel time and the cumulative 
%numbers of the vehicles departures and arrivals after we run the
%cell-transmission-model (CTM) from Daganzo in the function 
%CTM_network_load.m. 
%[act_travel_time,Arriving,Departing] =
%travel_time(x_veh,flow,paths,discretesplits_nodes,dt,source,num_demands).
%--------------------------------------------------------------------------
%inputs:
%x_veh: (timesteps + 1 x paths x number of cells)-matrix that shows the 
%maximum number of vehicles in cell i at time t for all paths.
%flow is a (timesteps x paths)-matrix, that shows the traffic flow leaving
%origin r at time t via path p between OD pair r-s.
%paths: a (1 x paths)-vector, that shows all existing paths.
%discretesplits_nodes: (2 x number of cells) gives us a correspondance 
%from each link to his beginning and ending node.
%dt: time interval for the simulation.
%source: source link.
%num_demands: Number of available time-steps for demand values.
%
% by Manuel Jakob
% 19 March 2012
%==========================================================================

T = size(x_veh,1)-1;
num_cells = size(x_veh,2);
n_paths = size(x_veh,3);
traffic_flow = zeros(T,num_cells,n_paths);

num_nodes = max(max(discretesplits_nodes));
%Initialize:
Departing = zeros(T,num_nodes,num_nodes);
Departing_paths = zeros(T,n_paths);
Arriving = zeros(T,num_nodes,num_nodes);
Arriving_paths = zeros(T,n_paths);
act_travel_time = zeros(num_demands,n_paths);

%Determine the traffic flow for all cells for all paths at all times 
%from the TOPL simulation tool output.
for t = 1:T %all times
    for p = 1:n_paths %all paths
        traffic_flow(t,1,p) = flow(t,p);
        cpath = cell2mat(paths(p));
            for I = 1:size(cpath,2)-2 %all nodes of path cpath
                for i = 1:num_cells-1 %all cells
                    if cpath(I) == discretesplits_nodes(1,i)
                        if discretesplits_nodes(2,i) == cpath(I+1)
                        
                            %Look for cell i+1 (the successor of cell i at
                            %path p)
                            for j = 1:num_cells %all cells
                                if cpath(I+1) == discretesplits_nodes(1,j)
                                    if discretesplits_nodes(2,j) == cpath(I+2)
                                    CellPlusOne = j;
                                    break
                                    end
                                end
                            end %all cells
                            traffic_flow(t,CellPlusOne,p) = traffic_flow(t,i,p) + x_veh(t,i,p) - x_veh(t+1,i,p);
                        break
                        end
                    end
                end %all cells
            end %all nodes of path cpath
    end %all paths 
end %all times

%Determine the cumulative numbers of the vehicles departures and arrivals 
%for each path at all times.
for p = 1:n_paths %all paths
    for t = 1:T %all times
        if t == 1
            Departing_paths(t,p) = traffic_flow(t,1,p)*dt;
            Arriving_paths(t,p) = traffic_flow(t,num_cells,p)*dt;
        else
            Departing_paths(t,p) = Departing_paths(t-1,p) + traffic_flow(t,1,p)*dt;
            Arriving_paths(t,p) = Arriving_paths(t-1,p) + traffic_flow(t,num_cells,p)*dt;
        end
    end %all times
end %all paths
  
%Determine the cumulative numbers of the vehicles departures and arrivals 
%for each cell at all times.
for t = 1:T-1 %all times
    for i = 1:num_cells %all cells
        I = discretesplits_nodes(1,i);
        J = discretesplits_nodes(2,i);
        if t == 1
            Departing(t+1,I,J) = sum(traffic_flow(t,i,:),3)*dt;
        else
            Departing(t+1,I,J) = Departing(t,I,J) + sum(traffic_flow(t,i,:),3)*dt;
        end
                
        traffic_flow_arr = 0;
        for j = 1:num_cells %all cells 
            K = discretesplits_nodes(2,j);
            %Look for output nodes of each cell that are equal to the input
            %node of cell i.
            if K == I
               traffic_flow_arr = traffic_flow_arr + sum(traffic_flow(t,j,:),3)*dt;
               if t == 1
                   Arriving(t+1,I,J) = traffic_flow_arr;
               else
                   Arriving(t+1,I,J) = Arriving(t,I,J) + traffic_flow_arr;
               end    
            end
        end %all cells
    end %all cells
end %all times
%the cumulative arriving counts for the source link.
for t = 1:T %all times
    for i = 1:num_cells %all cells
        if i == source
            I = discretesplits_nodes(1,i);
            J = discretesplits_nodes(2,i);
            if t == 1
                Arriving(t,I,J) = sum(flow(t,:),2)*dt;
            else
                Arriving(t,I,J) = Arriving(t-1,I,J) + sum(flow(t,:),2)*dt;
            end
            break
        end
    end %all cells
end %all times

%Determine actual travel time, by summing up the discrete vehicle values.
for t = 1:num_demands %all times
    for p = 1:n_paths %all paths
    
    if t == 1
        denominator = Departing_paths(t,p);
        Arr = Arriving_paths(Arriving_paths(:,p) <= Departing_paths(t,p) & Arriving_paths(:,p) > 0,p);
        if isempty(Arr)
            act_travel_time(t,p) = 0;
        else
        nominator = ...
        sum((find(Arriving_paths(:,p) <= Departing_paths(t,p) & Arriving_paths(:,p) > 0,size(Arr,1)) - t)...
        .*diff([0; Arr]));
    
        last_entry_Arr = Arriving_paths(find(Arriving_paths(:,p) <= Departing_paths(t,p) & Arriving_paths(:,p) > 0,1,'last'),p);
        if last_entry_Arr < Departing_paths(t,p) 
            rest = (Departing_paths(t,p) - last_entry_Arr)*(find(Arriving_paths(:,p) <= Departing_paths(t,p) & Arriving_paths(:,p) > 0,1,'last')+1 - t);
            nominator = nominator + rest;
        end
        act_travel_time(t,p) = nominator/denominator*dt;
        end
    else   
        denominator = Departing_paths(t,p) - Departing_paths(t-1,p);
        Arr = Arriving_paths(Arriving_paths(:,p) <= Departing_paths(t,p) & Arriving_paths(:,p) > Departing_paths(t-1,p),p);
        if isempty(Arr)
            act_travel_time(t,p) = 0;
        else
        nominator = ...
        sum((find(Arriving_paths(:,p) <= Departing_paths(t,p) & Arriving_paths(:,p) > Departing_paths(t-1,p),size(Arr,1)) - t)...
        .*diff([Departing_paths(t-1,p); Arr]));
        
        last_entry_Arr = Arriving_paths(find(Arriving_paths(:,p) <= Departing_paths(t,p) & Arriving_paths(:,p) > Departing_paths(t-1,p),1,'last'),p);
        if last_entry_Arr < Departing_paths(t,p) 
            rest = (Departing_paths(t,p) - last_entry_Arr)*(find(Arriving_paths(:,p) <= Departing_paths(t,p) & Arriving_paths(:,p) > Departing_paths(t-1,p),1,'last')+1 - t);
            nominator = nominator + rest;
        end 
        act_travel_time(t,p) = nominator/denominator*dt;
        end
    end
    
    end %all paths
end %all times

% %Determine actual travel time, using the trapezoidal rule.
% for t = 1:num_demands %all times
%     for p = 1:n_paths %all paths
%         %Determine A^(-1)(D(t)) and A^(-1)(D(t-1)) for time t and for path p.
%         [C_t] = Inverse_calc(Arriving_paths,Departing_paths,t,p);
%         if t > 1 
%             [C_tm1] = Inverse_calcm1(Arriving_paths,Departing_paths,t,p);
%             %Numerical integration, using trapezoidal rule, t > 1.
%             act_travel_time(t,p) = (C_t + C_tm1 - 2*t + 1)*dt/2;
%         else
%             %Numerical integration, using trapezoidal rule, t = 1.
%             act_travel_time(t,p) = (C_t - 2*t + 1)*dt/2;
%         end
%     end %all paths
% end %all times

end

% -------------------------------------------------------------------------
% function [C_t] = Inverse_calc(Arriving_paths,Departing_paths,t,p)
%     %Interpolation using the cumulative number of arrivals and departures
%     %to compute the cost matrix of A^(-1)(D(t)) for time t and for path p,
%     %where D is the cumulative number of departures and A is the
%     %cumulative number of arrivals for time t and for path p.
%     
%     Dep = Departing_paths(t,p);
%     interpolA = find(Arriving_paths(:,p) <= Dep,1,'last');
%     if isempty(interpolA)
%         interpolA = 1;
%     end
%     vehA = Arriving_paths(interpolA,p);
%     interpolB = find(Arriving_paths(:,p) >= Dep,1,'first');
%     if isempty(interpolB)
%         interpolB = size(Arriving_paths,1);
%     end
%     vehB = Arriving_paths(interpolB,p);
%     delta = (Dep - vehA)/(vehB - vehA);
%     if vehB - vehA == 0
%         C_t = interpolB;
%     else
%         C_t = interpolA + delta*(interpolB - interpolA);
%     end
%      
% end
% 
% % -------------------------------------------------------------------------
% function [C_tm1] = Inverse_calcm1(Arriving_paths,Departing_paths,t,p)
%     %Interpolation using the cumulative number of arrivals and departures
%     %to compute the cost matrix of A^(-1)(D(t-1)) for time t and for path p,
%     %where D is the cumulative number of departures and A is the
%     %cumulative number of arrivals for time t and for path p.
%        
%     Dep = Departing_paths(t-1,p);
%     interpolA = find(Arriving_paths(:,p) <= Dep,1,'last');
%     if isempty(interpolA)
%         interpolA = 1;
%     end
%     vehA = Arriving_paths(interpolA,p);
%     interpolB = find(Arriving_paths(:,p) >= Dep,1,'first');
%     if isempty(interpolB)
%         interpolB = size(Arriving_paths,1);
%     end
%     vehB = Arriving_paths(interpolB,p);
%     delta = (Dep - vehA)/(vehB - vehA);
%     if vehB - vehA == 0
%         C_tm1 = interpolB;
%     else
%         C_tm1 = interpolA + delta*(interpolB - interpolA); 
%     end
%     
% end


