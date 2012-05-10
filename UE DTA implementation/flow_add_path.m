function [flow,paths] = flow_add_path(npaths_f,flow,iteration,t,demand,source,paths,shortest_paths,xml_file_Java,diff_length)

%When the auxiliary path is newly generated (shortest_paths(t) is 
%not existing in the matrix paths), flow_add_path is adding a new path
%entry into the vector paths and the corresponding values into the flow
%matrix.
%[flow,paths] =
%flow_add_path(npaths_f,flow,iteration,t,...
%demand,source,paths,shortest_paths,xml_file_Java,diff_length)
%--------------------------------------------------------------------------
%inputs:
%npaths_f: number of paths.
%flow is a (timesteps x paths)-matrix, that shows the traffic flow leaving
%origin r at time t via path p between OD pair r-s.
%iteration: current step of iterations
%t: current time step t.
%demand: compliant demand-vector plus noncompliant demand-vector for
%all timesteps.
%source: source cell.
%paths is a (1 x paths)-vector, that shows all existing paths.
%shortest_paths is a (1 x timesteps)-vector, that shows the shortest paths
%of the given network, created by the Dijkstra algorithm.
%xml_file_Java: name of the (.xml) network file for Java.
%diff_length = 1, when we consider paths with different length and = 0,
%when paths has the same length (number of nodes of the paths are equal).
%
% by Manuel Jakob
% 06 May 2012
%==========================================================================


for i = 1:npaths_f %for all paths
   flow(t,i) = iteration/(iteration + 1)*flow(t,i);
end %for all paths
flow(t,npaths_f + 1) = 1/(iteration + 1)*demand(t,source); 

if diff_length == 1 %if we consider paths with different length
    
for pp = 1:size(paths,2)
   if size(cell2mat(shortest_paths(t)),2) == size(cell2mat(paths(pp)),2)
       break
   else
       if pp == size(paths,2)
           paths(pp + 1) = shortest_paths(t);
           %Changes the (.xml)-file due to the additional number of paths:
           myXMLEditor(xml_file_Java)
       end
   end
end

else %paths has the same length (number of nodes of the paths are equal) 
    
for pp = 1:size(paths,2)
   if cell2mat(shortest_paths(t)) == cell2mat(paths(pp))
       break
   else
       if pp == size(paths,2)
           paths(pp + 1) = shortest_paths(t);
           %Changes the (.xml)-file due to the additional number of paths:
           myXMLEditor(xml_file_Java)
       end
   end
end

end %if we consider paths with different length


