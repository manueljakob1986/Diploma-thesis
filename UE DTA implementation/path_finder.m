function [corr_paths] = path_finder(I,J,paths)

%path_finder determines the corresponding paths, in which the nodes I and J
%exist.
%[corr_paths] = path_finder(I,J,paths)
%--------------------------------------------------------------------------
%inputs:
%I: starting node of a link.
%J: ending node of a link.
%paths: a (1 x paths)-vector, that shows all existing paths.
%
% by Manuel Jakob
% 30 April 2012
%==========================================================================

npaths = size(paths,2);
i = 1;
corr_paths = 0;

for p = 1:npaths %for all paths
one_element = 0;
two_element = 0;
path = paths{p};

    for K = 1:size(path,2) %for all nodes in path
        if path(K) == I
            one_element = 1;
        end
        if path(K) == J
            two_element = 1;
        end
    end %for all nodes in path
    if one_element == 1 && two_element == 1
        corr_paths(i) = p;
        i = i + 1;
    end
    
end %for all paths


