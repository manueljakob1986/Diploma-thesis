function [jSplitRatioMatrix] = java_three_translate(paths_TR)

%Here we traslate the three-dimensional matrix paths_TR into Java input
%with data type Double.
%[jSplitRatioMatrix] = Java_three_translate(paths_TR)
%--------------------------------------------------------------------------
%inputs:
%paths_TR: Turning ratio matrix for all paths for all nodes (with an 
%input and output link) at all times.
%
% by Manuel Jakob
% 18 April 2012
%==========================================================================

%Translate paths_TR for Java input.
[~,numLinks,numPaths] = size(paths_TR);
jSplitRatioMatrix = javaArray('java.lang.Double[][]', numLinks);

for iLinksOne = 1:numLinks
   jLinksTwo = javaArray('java.lang.Double[]', numLinks);
       for iLinksTwo = 1:numLinks
           jPath = javaArray('java.lang.Double', numPaths);
               for iPath = 1:numPaths
                   jPath(iPath) = java.lang.Double(paths_TR(iLinksOne,iLinksTwo,iPath));
               end
               jLinksTwo(iLinksTwo) = jPath;
       end
       jSplitRatioMatrix(iLinksOne) = jLinksTwo;
end
