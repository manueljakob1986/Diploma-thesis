function [new_flow,new_paths] = flow_update(iteration,flow,paths,shortest_paths,demand,source,xml_file_Java,num_demands)

%flow_update determines the updated flow explained in the paper from Tong
%and Wong, 1999.
%[new_flow,new_paths] = flow_update(iteration,flow,paths,...
%shortest_paths,demand,source,xml_file_Java,num_demands).
%--------------------------------------------------------------------------
%inputs:
%iteration: number of interations.
%flow is a (timesteps x paths)-matrix, that shows the traffic flow leaving
%origin r at time t via path p between OD pair r-s.
%paths is a (1 x paths)-vector, that shows all existing paths.
%shortest_paths is a (1 x timesteps)-vector, that shows the shortest paths
%of the given network, created by the Dijkstra algorithm.
%demand: compliant demand-vector plus noncompliant demand-vector for
%all timesteps.
%source: source cell.
%xml_file_Java: name of the (.xml) network file for Java.
%num_demands: Number of available time-steps for demand values.
%
% by Manuel Jakob
% 12 March 2012
%==========================================================================

[~,npaths_f] = size(flow);
npaths_p = size(paths,2);

%cases that are not possible:
if npaths_f ~= npaths_p 
    display('This is not a valid input.')
    return
end

%Update the flow-matrix, based on the flow-updating scheme in the Tong,
%Wong paper,1999.
for t = 1:num_demands %for all timesteps
   for p = 1:npaths_f %for all paths
       
       %Check whether shortest_paths and paths have the same length (same
       %number of nodes in the paths entries)
       if size(cell2mat(shortest_paths(t)),2) == size(cell2mat(paths(p)),2) %check for length
   
       if cell2mat(shortest_paths(t)) == cell2mat(paths(p))
       %if the auxiliary path is an old path (shortest_paths(t) is 
       %already existing in the matrix paths).
           save_path = p;
           for i = 1:npaths_f %for all paths
               if i == save_path
                   flow(t,i) = iteration/(iteration + 1)*flow(t,i) + 1/(iteration + 1)*demand(t,source);  
               else
                   flow(t,i) = iteration/(iteration + 1)*flow(t,i);
               end
           end %for all paths
           break
       else
       %The auxiliary path is newly generated (shortest_paths(t) is 
       %not existing in the matrix paths).
           if p == npaths_f
               [flow,paths] = flow_add_path(npaths_f,flow,iteration,t,demand,source,paths,shortest_paths,xml_file_Java,0);
           end
       end
       else
           if p == npaths_f
               [flow,paths] = flow_add_path(npaths_f,flow,iteration,t,demand,source,paths,shortest_paths,xml_file_Java,1);
           end
       end %check for length
   end %for all paths
end %for all timesteps
%Update the matrices flow and paths.
new_flow = flow;
new_paths = paths;

end
