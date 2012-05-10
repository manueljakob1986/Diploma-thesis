function [flow_result,accuracy,accuracy1,accuracy2,flow,iteration,...
    act_tt_result,num_demands] = uedta(xml_file,epsilon,iteration_max,xml_file_Java)

%Here we run the UE DTA iterative algorithm. The algorithm stopps, when
%the maximum number of iterations iteration_max is reached or when the
%accuracy of the convergence criteria is sufficient.
%[flow,iteration] = uedta(xml_file,epsilon,iteration_max,xml_file_Java).
%--------------------------------------------------------------------------
%inputs:
%xml_file: name of the (.xml) file.
%epsilon: stopping criteria.
%iteration_max: maximum number of iterations.
%xml_file_Java: name of the (.xml) network file for Java.
%
% by Manuel Jakob
% 12 March 2012
%==========================================================================

%initialize the multiplicator of the time horizon T.
T_multiplicator = 1;

%define input data with the help of the (.xml) read function topl_myxmlread_nosplit.m.
[T, num_cells, inc_matrix, ~, ~, ~, compliant_demand, ...
    noncompliant_demand, ~, ~, source, sink, dt] ...
    = topl_myxmlread_nosplit(xml_file, T_multiplicator);

%Determine the initial_flow, vector of dimension (timesteps x 1).
total_demand = sum(compliant_demand,2) + sum(noncompliant_demand,2);
demand = compliant_demand + noncompliant_demand;
initial_flow = total_demand;
%Number of available time-steps for demand values.
num_demands = find(initial_flow(:) <= 0,1,'first')-1;

%discretesplits_nodes gives us a correspondance from each link to his
%beginning and ending node.
[~, discretesplits_nodes, ~, freeflowspeed, length] = splitlink(xml_file);
num_nodes = max(max(discretesplits_nodes));

%Determine SID and FID.
SID = discretesplits_nodes(1,source);
FID = discretesplits_nodes(2,sink);

%Determine the initial_paths, vector of dimension (timesteps x 1):
Arriving = zeros(T,num_nodes,num_nodes);
Departing = zeros(T,num_nodes,num_nodes);
%Run the single-origin-single-destination shortest paths Dijkstra
%algorithm.
initial_case = 1;
flowratio = ones(num_demands,1);
paths = num2cell(nan(1));
travel_times = zeros(T,num_cells);
 for i = 1:num_cells %all cells
    travel_times(:,i) = length(i)/freeflowspeed(i)*60;
 end %all cells
[~,shortest_paths] = td_dijkstra(inc_matrix,Arriving,Departing,SID,FID,dt,...
    discretesplits_nodes,initial_case,travel_times,num_demands,flowratio,paths);
initial_paths = shortest_paths(1);
initial_case = 0;

%initialize the matrices paths and flow:
%initial_paths is a (timesteps x 1)-matrix, that shows the intial paths for 
%all timesteps t.
paths = initial_paths;
%initial_flow is a (timesteps x 1)-matrix, that shows the initial flow 
%leaving origin r via the intial path initial_paths between the OD pair
%r-s for all timesteps t.
flow = initial_flow;

%UE DTA iterative algorithm, maximum number of iterations is iteration_max.
for iteration = 1:iteration_max

    %Load the network due to the cell-transmission model in TOPL.
    [x_veh] = CTM_network_load(dt,T,flow,paths,discretesplits_nodes);
    
    %Compute the actual travel time and the cumulative numbers of 
    %the vehicles departures and arrivals.
    [act_travel_time,Arriving,Departing] = travel_time(x_veh,flow,paths,...
        discretesplits_nodes,dt,source,num_demands);
    
    %Determine the convergence criteria, i.e., the value of the 
    %duality gap of the objective function from the predictive dynamic user 
    %equilibrium minimisation problem.
    %And we calculate the shortest paths of the given network for all times, 
    %created by the time-dependent single-origin-single-destination Dijkstra 
    %algorithm.
    [G,G1,G2,shortest_paths] = convergence_criteria(flow,act_travel_time,inc_matrix,...
        Arriving,Departing,SID,FID,dt,discretesplits_nodes,initial_case,...
        travel_times,num_demands,flowratio,paths);
    
    
    %Save the resulting flow for the first two time steps and the accuracy 
    %at each UE DTA iteration.
    for p = 1:size(paths,2)

        for t = 1:num_demands
            flow_result(t,iteration,p) = flow(t,p);
            act_tt_result(t,iteration,p) = act_travel_time(t,p);
        end
    end
        accuracy(iteration) = G; %total accuracy for all times
        accuracy1(iteration) = G1; %accuracy at time step 1
        accuracy2(iteration) = G2; %accuracy at time step 2
    
    
%     %if stopping criteria is sufficient, then stop.
%     if G < epsilon
%         display('Accuracy of the convergence criteria is sufficient. UE DTA algorithm is stopped.')
%         break
%     end
    
    
    %Update the network flow.
    [flow,paths] = flow_update(iteration,flow,paths,shortest_paths,demand,source,...
        xml_file_Java,num_demands);
    
    
    %Determines the flow ratio from one path p at time t relative to the flow on all paths.
    [flowratio] = flow_ratio(num_demands,flow);
end
end
