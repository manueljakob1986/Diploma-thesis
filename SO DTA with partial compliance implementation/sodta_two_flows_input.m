% Here we run the function sodta_two_flows.m with the input created by the
% function [] = topl_myxmlread(fileName, T_multiplicator).
% sodta_two_flows.m solves the single destination SO-DTA LP.
%
% by Manuel Jakob
% 06 February 2012
%==========================================================================


clc; clear all; close all;

%initialize the multiplicator of the time horizon T
t_mult_start = 1;

%define input data for the function sodta_two_flows.m with the function topl_myxmlread.m
[T, num_cells, inc_matrix, q_max, deltas, N_max, compliant_demand, noncompliant_demand, ...
    turning_ratio, compliant_initial_vehicles, noncompliant_initial_vehicles,~,sink,~] = ...
    topl_myxmlread('mySplitNetwork_order.xml', t_mult_start);

%run function sodta_two_flows.m
tic
[simulation, T_multiplicator] = sodta_two_flows(T, num_cells, sink, inc_matrix, q_max, deltas, N_max, ...
    compliant_demand, noncompliant_demand, turning_ratio, compliant_initial_vehicles, noncompliant_initial_vehicles, t_mult_start);
toc

%repeat simulation with a doubled time-horizon T, when the total inflow and the total outflow of
%the network are not equal
while simulation ~= 0
    if simulation == 1   
        clc; close all;
        disp(['The new time horizon is ',num2str(T_multiplicator),'-times the old time horizon.'])
        
        [T, num_cells, inc_matrix, q_max, deltas, N_max, compliant_demand, noncompliant_demand, ...
        turning_ratio, compliant_initial_vehicles, noncompliant_initial_vehicles,~,sink,~] = ...
        topl_myxmlread('mySplitNetwork_order.xml', T_multiplicator);
        
        [simulation, T_multiplicator] = sodta_two_flows(T, num_cells, sink, inc_matrix, q_max, deltas, N_max, ...
           compliant_demand, noncompliant_demand, turning_ratio, compliant_initial_vehicles, noncompliant_initial_vehicles, T_multiplicator);
    end
end
