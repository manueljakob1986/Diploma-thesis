function [sp_travel_time,paths] = td_dijkstra(incidence_matrix,Arriving,Departing,...
    SID,FID,dt,discretesplits_nodes,initial_case,travel_times,num_demands,flowratio,all_paths)

%time-dependent single-origin-single-destination DIJKSTRA Calculate Minimum 
%Costs and Paths using time-dependent Dijkstra's Algorithm
%[sp_travel_time,paths] = td_dijkstra(incidence_matrix,Arriving,Departing,...
%    SID,FID,dt,discretesplits_nodes,initial_case,travel_times,...
%    num_demands,flowratio,all_paths).
%--------------------------------------------------------------------------
%inputs:
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

[n,nc] = size(incidence_matrix);
[at,ac,acc] = size(Arriving);
[Dt,Dc,Dcc] = size(Departing);
[E] = incidence(discretesplits_nodes);

%cases that are not possible:
if at ~= Dt 
    display('This is not a valid input.')
    return
end
if ac ~= Dc 
    display('This is not a valid input.')
    return
end
if n ~= nc 
    display('This is not a valid input.')
    return
end
if ac ~= acc 
    display('This is not a valid input.')
    return
end
if Dc ~= Dcc 
    display('This is not a valid input.')
    return
end

% Find the Minimum Costs and Paths using Dijkstra's Algorithm
% Initializations
min_cost = Inf(1,ac);
sp_travel_time = zeros(1,num_demands);
path = num2cell(nan(1,ac));
paths = num2cell(nan(1,num_demands));
gamma = zeros(num_demands,ac);
for t = 1:num_demands
    for i = 1:ac
        if i == SID
           gamma(t,i) = t*dt;
        else
           gamma(t,i) = Inf;   
        end
    end
end
table = zeros(1,ac);

%for all timesteps
for t = 1:num_demands
I = SID;
table(I) = 1;
path(I) = {I};
    
while any(table)
nids = find(E(:,1) == I);
% Calculate the Costs to the Neighbor Points and Record Paths
for kk = 1:length(nids)
    vehtime = gamma(t,I)/dt;
    J = E(nids(kk),2);
    
    %If we have more than one output cells (splitting point), compute the
    %flow ratio from one path p at time t relative 
    %to the flow on all paths.
    FR = 0;
    if length(nids) > 1 && size(all_paths,2) > 1
        [corr_paths] = path_finder(I,J,all_paths);
        if corr_paths == 0
            FR = 0;
        else
            for cp = 1:size(corr_paths,2)
                FR = FR + flowratio(t,corr_paths(cp));
            end
        end
    end

    if initial_case == 1
        [C] = processInputs_travel_time(travel_times,discretesplits_nodes,vehtime,t,I,J);
         C = C*dt;
    else %if initial_case == 0
        [C] = processInputs(Arriving,Departing,vehtime,I,J,FR,travel_times,discretesplits_nodes,t);
         C = C*dt;
    end
    
    if C < gamma(t,J)
        gamma(t,J) = C;
        path{J} = [path{I} J];
        if J ~= FID
            table(J) = 1;
        end
    end
end

if J ~= FID
    table(I) = 0;
else 
    table = zeros(1,ac);
end

nodes = find(table);
N = find(gamma(t,nodes) == min(gamma(t,nodes)),1);

if isempty(N)
% Settle the last Minimum Value for each timestep
min_cost(J) = gamma(t,I);     
else
% Settle the Minimum Value
I = nodes(N);
min_cost(I) = gamma(t,I);
end

end
% Store actual travel time and Paths
sp_travel_time(t) = min_cost(FID) - t*dt;
paths(t) = path(FID);

end    
end
% -------------------------------------------------------------------------
function [C] = processInputs(Arriving,Departing,vehtime,I,J,FR,travel_times,discretesplits_nodes,t)
    %interpolation using the cumulative number of arrivals and departures
    %to compute the cost matrix of D^(-1)(A(t)) for all links and all times
    %t, where D is the cumulative number of departures and A is the
    %cumulative number of arrivals for all links and all times

    beta = (vehtime - floor(vehtime))/(ceil(vehtime) - floor(vehtime));
    if ceil(vehtime) - floor(vehtime) == 0
        beta = 0;
    end
    if vehtime > size(Arriving,1)
       vehtime = size(Arriving,1);
    end
    Arr = Arriving(floor(vehtime),I,J) + beta*(Arriving(ceil(vehtime),I,J) ...
        - Arriving(floor(vehtime),I,J));
    if FR ~= 0
        Arr = Arr * FR;
    end
       
    interpolA = find(Departing(:,I,J) <= Arr,1,'last');
    if isempty(interpolA)
        interpolA = 1;
    end
    vehA = Departing(interpolA,I,J);
    
    interpolB = find(Departing(:,I,J) >= Arr,1,'first');
    if isempty(interpolB)
        interpolB = size(Departing,1);
    end
    vehB = Departing(interpolB,I,J);
    delta = (Arr - vehA)/(vehB - vehA);
    if vehB - vehA == 0
        C = interpolB;
    else
        C = interpolA + delta*(interpolB - interpolA); 
    end
    
    if Departing(:,I,J) == zeros(size(Departing,1),1)
        [C] = processInputs_travel_time(travel_times,discretesplits_nodes,vehtime,t,I,J);
    end
end
% -------------------------------------------------------------------------
function [TTcum] = processInputs_travel_time(travel_times,discretesplits_nodes,TTcum,t,I,J)
%Determine the cumulative travel time for link i (nodes I and J). This
%case is only used for the initial flow time-dependent Dijkstra algorithm.

for i = 1:size(discretesplits_nodes,2) %for all links
    if I == discretesplits_nodes(1,i)
        if J == discretesplits_nodes(2,i)
            link = i;
            TTcum = TTcum + travel_times(t,link);
            break
        end
    end
end %for all links
end
% -------------------------------------------------------------------------
% Output Edge List
function [E] = incidence(discretesplits_nodes)  
    E = discretesplits_nodes';
end

       