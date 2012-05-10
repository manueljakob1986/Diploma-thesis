function [G,G1,G2,shortest_paths] = convergence_criteria(flow,act_travel_time,...
    incidence_matrix,Arriving,Departing,SID,FID,dt,discretesplits_nodes,...
    initial_case,travel_times,num_demands,flowratio,all_paths)

%cost_function determines the convergence criteria, i.e. the value of the 
%duality gap of the objective function from the predictive dynamic user 
%equilibrium minimisation problem.
%And we calculate the shortest paths of the given network for all times, 
%created by the time-dependent single-origin-single-destination Dijkstra 
%algorithm.
%[G,shortest_paths] = convergence_criteria(flow,act_travel_time,...
%    incidence_matrix,Arriving,Departing,SID,FID,dt,discretesplits_nodes,...
%    initial_case,travel_times,num_demands,flowratio,all_paths).
%--------------------------------------------------------------------------
%inputs:
%flow is a (timesteps x paths)-matrix, that shows the traffic flow leaving
%origin r at time t via path p between OD pair r-s.
%act_travel_time is a (timesteps x paths)-matrix, that shows the actual
%travel time on path p at time t between OD pair r-s.
%incidence_matrix is the incidence matrix, we get from the network editor
%(.xml) file.
%Arriving is a (timesteps x vertices x vertices)-matrix that shows 
%cumulative numbers of the vehicles arrivals for all timesteps and for all links.
%Departing is a (timesteps x vertices x vertices)-matrix that shows 
%cumulative numbers of the vehicles departures for all timesteps and for all links.
%SID is a the starting point (we only consider single origin problems).
%FID is the finish point (we only consider single destination problems).
%dt: time interval for the simulation.
%discretesplits_nodes: (2 x number of cells) gives us a correspondance 
%from each link to his beginning and ending node.
%initial_case: = 1 for the initial flow case, otherwise equal to 0.
%travel_times: (timesteps x cells)-matrix that shows the travel times for
%the vehicles at each cell at all times.
%num_demands: Number of available time-steps for demand values.
%flowratio: (num_demands x paths) matrix showing the flow ratio from one
%path p at time t relative to the flow on all paths.
%all_paths is a (1 x paths)-vector, that shows all existing paths.
%
% by Manuel Jakob
% 08 March 2012
%==========================================================================

[~,paths] = size(flow);

%Run the single-origin-single-destination shortest paths Dijkstra
%algorithm:
%Output: 
%sp_travel_time is a (1 x timesteps)-vector, that shows the minimum
%travel time at time t between OD pair r-s.
%shortest_paths is a (1 x timesteps)-vector, that shows the shortest paths
%of the given network, created by the Dijkstra algorithm.
[~,shortest_paths] = td_dijkstra(incidence_matrix,Arriving,...
    Departing,SID,FID,dt,discretesplits_nodes,initial_case,travel_times,...
    num_demands,flowratio,all_paths);

%Determine the demand-vector of dimension timesteps by using the
%flow-matrix for all times.
demand = zeros(1,num_demands);
sp_travel_time = zeros(1,num_demands);
for t = 1:num_demands
    demand(t) = sum(flow(t,:),2);
    sp_travel_time(t) = min(act_travel_time(t,:));
end

%Determine the convergence criteria.
G = sum(sum(flow(1:num_demands,:).*abs(act_travel_time - sp_travel_time'*ones(1,paths))))/(demand*sp_travel_time')*100;
%convergence criteria for time step 1
G1 = sum(sum(flow(1,:).*abs(act_travel_time(1,:) - sp_travel_time(1)'*ones(1,paths))))/(demand(1)*sp_travel_time(1)')*100;
%convergence criteria for time step 2
G2 = sum(sum(flow(num_demands,:).*abs(act_travel_time(num_demands,:) - sp_travel_time(num_demands)'*ones(1,paths))))/(demand(num_demands)*sp_travel_time(num_demands)')*100;


end