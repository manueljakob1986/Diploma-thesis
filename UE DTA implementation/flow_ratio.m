function [flowratio] = flow_ratio(num_demands,flow)

%flow_ratio determines the flow ratio from one path p at time t relative 
%to the flow at time t on all paths.
%[flowratio] = flow_ratio(num_demands,flow)
%--------------------------------------------------------------------------
%inputs:
%num_demands: Number of available time-steps for demand values.
%flow is a (timesteps x paths)-matrix, that shows the traffic flow leaving
%origin r at time t via path p between OD pair r-s.
%
% by Manuel Jakob
% 30 April 2012
%==========================================================================

n_paths = size(flow,2);
flowratio = zeros(num_demands,n_paths);
for t = 1:num_demands
    for p = 1:n_paths
        flowratio(t,p) = flow(t,p)/sum(flow(t,:),2);
    end
end