function [paths_TR] = paths_turning_ratios(paths,discretesplits_nodes)

%We determine the corresponding turning ratio matrix paths_TR for all paths
%for all nodes (with an input and output link) at all times.
%[paths_TR] = paths_turning_ratios(paths,discretesplits_nodes).
%--------------------------------------------------------------------------
%inputs:
%paths is a (1 x paths)-vector, that shows all existing paths.
%discretesplits_nodes: (2 x number of cells) gives us a correspondance 
%from each link to his beginning and ending node.
%
% by Manuel Jakob
% 17 April 2012
%==========================================================================

num_paths = size(paths,2);
num_links = size(discretesplits_nodes,2);
%Initialize:
paths_TR = zeros(num_links,num_links,num_paths);

for p = 1:num_paths %for all paths
path = cell2mat(paths(p));
    for i = 1:size(path,2) %for all nodes of path
        for j = 1:num_links %for all links
            %determination of the input link
            if discretesplits_nodes(2,j) == path(i) %link input
                if discretesplits_nodes(1,j) == path(i-1)
                    input_link = j;
                end
                for k = 1:num_links %for all links
                    %determination of the output link
                    if discretesplits_nodes(1,k) == path(i) %link output
                        if discretesplits_nodes(2,k) == path(i+1)
                            output_link = k;
                            paths_TR(input_link,output_link,p) = 1;
                        end
                    end %link output
                end %for all links
            end %link input
        end %for all links
    end %for all nodes of path
end %for all paths



