%Here you can run the example network in the Ziliaskopoulus paper (Figure
%3) or the parallel flow network with an alternate to the main line in
%different versions.
%First you choose between the following networks:
%A: Example Ziliaskopoulus network
%B: Parallel flow network with two destinations
%C: Parallel flow network with one destinations.
%
%After that you choose whether you want to have only compliant flow demand 
%or compliant and non-compliant flow demand in the network.
%The total demand in the network stays equal in both variants
%,i.e., second you choose the input demand:
%A: Only compliant flow demand
%B: Compliant and non-compliant flow demand.
%
% by Manuel Jakob
% 29 November 2011


% define general input data for the function sodta_two_flows.m:
num_cells = 10;
T = 20;
deltas = ones(T,num_cells);

N_max = zeros(T,num_cells);
for i = 1:T
    N_max(i,:) = [10000,20,10,10,10,10,10,10,20,10000];
end    

n_cars_init_comp = [0,0,0,0,0,0,0,0,0,0];
n_cars_init_noncomp = [0,0,0,0,0,0,0,0,0,0];

%==========================================================================
%Choose one network:
network = 0;
destination = 0;
reply = input('Please choose one of the following networks: \n Example Ziliaskopoulus network [A] \n Parallel flow network with two destinations  [B] \n Parallel flow network with one destination  [C] \n (default value [A]): ', 's');
if isempty(reply)
    reply = 'A';
end
if reply == 'A'
    network = 1;
elseif reply == 'B'
    network = 2;
elseif reply == 'C'
    network = 3;
end
%==========================================================================
if network == 1
    %Example Ziliaskopoulus network.
       
inc_matrix = ...
   [0,1,0,0,0,0,0,0,0,0 ;
    0,0,1,0,1,0,0,0,0,0 ;
    0,0,0,1,0,1,0,0,0,0 ;
    0,0,0,0,0,0,0,0,1,0 ;
    0,0,0,0,0,0,1,0,0,0 ;
    0,0,0,0,0,0,1,0,0,0 ;
    0,0,0,0,0,0,0,1,0,0 ;
    0,0,0,0,0,0,0,0,1,0 ;
    0,0,0,0,0,0,0,0,0,1 ;
    0,0,0,0,0,0,0,0,0,0];

q_max = ...
   [12,12,6,6,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000;
    12,12,6,0,6,6,6,6,12,10000;
    12,12,6,0,6,6,6,6,12,10000;
    12,12,6,3,6,6,6,6,12,10000;
    12,12,6,3,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000;
    12,12,6,6,6,6,6,6,12,10000];

turning_ratios = [0,0.5,0.5,0,0,0,0,0,0,0];
destination = 1;
end
%==========================================================================
if network == 2
    %Parallel flow network with two destinations.
   
inc_matrix = ...
   [0,1,0,0,0,0,0,0,0,0 ;
    0,0,1,1,0,0,0,0,0,0 ;
    0,0,0,0,0,0,0,0,1,0 ;
    0,0,0,0,0,0,0,1,0,0 ;
    0,0,0,0,0,0,0,0,1,0 ;
    0,0,0,0,0,0,1,0,0,0 ;
    0,0,0,0,0,0,0,1,0,0 ;
    0,0,0,0,1,0,0,0,0,1 ;
    0,0,0,0,0,0,0,0,0,0 ;
    0,0,0,0,0,0,0,0,0,0];

q_max = zeros(T,num_cells);
for i = 1:T
    q_max(i,:) = [12,12,6,6,6,6,6,6,12,10000];
end 

turning_ratios = [0,0.5,0,0,0,0,0,0.5,0,0];
destination = 2;
end
%==========================================================================
if network == 3
    %Parallel flow network with one destination.
  
inc_matrix = ...
   [0,1,0,0,0,0,0,0,0,0 ;
    0,0,1,1,0,0,0,0,0,0 ;
    0,0,0,0,0,0,0,0,1,0 ;
    0,0,0,0,0,0,0,1,0,0 ;
    0,0,0,0,0,0,0,0,1,0 ;
    0,0,0,0,0,0,1,0,0,0 ;
    0,0,0,0,0,0,0,1,0,0 ;
    0,0,0,0,1,0,0,0,0,0 ;
    0,0,0,0,0,0,0,0,0,1 ;
    0,0,0,0,0,0,0,0,0,0];

q_max = zeros(T,num_cells);
for i = 1:T
    q_max(i,:) = [1000,6,5,5,5,5,5,5,6,10000];
end 

turning_ratios = [0,0.8,0,0,0,0,0,0,0,0];
destination = 1;
end
%==========================================================================
%Choose one demand:
flow_demand = 0;
reply = input('Please choose one of the following demands: \n Only compliant flow demand [A] \n Compliant and non-compliant flow demand  [B] \n (default value [A]): ', 's');
if isempty(reply)
    reply = 'A';
end
if reply == 'A'
    flow_demand = 1;
elseif reply == 'B'
    flow_demand = 2;
end
%==========================================================================
if flow_demand == 1
    %Only compliant flow demand.
    
compliant_demand = zeros(size(N_max));
compliant_demand(1,1) = 52;
compliant_demand(2,1) = 52;
compliant_demand(3,1) = 52;

noncompliant_demand = zeros(size(N_max));

end      
%==========================================================================
if flow_demand == 2
    %Compliant and non-compliant flow demand.
    
compliant_demand = zeros(size(N_max));
compliant_demand(1,1) = 26;
compliant_demand(2,1) = 26;
compliant_demand(3,1) = 26;

noncompliant_demand = zeros(size(N_max));
noncompliant_demand(1,1) = 26;
noncompliant_demand(2,1) = 26;
noncompliant_demand(3,1) = 26;

end 
%==========================================================================    
 
sodta_two_flows(T, num_cells, destination, inc_matrix, q_max, deltas, N_max, compliant_demand, noncompliant_demand, turning_ratios, n_cars_init_comp, n_cars_init_noncomp);
