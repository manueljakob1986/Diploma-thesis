function [jSourcedemand] = java_two_translate(sourcedemand)

%Here we traslate the two-dimensional matrix sourcedemand into Java input
%with data type Double.
%[jSourcedemand] = java_two_translate(sourcedemand)
%--------------------------------------------------------------------------
%inputs:
%sourcedemand: demand of the source cell for all paths.
%
% by Manuel Jakob
% 18 April 2012
%==========================================================================

%Translate sourcedemand for Java input.
[numTimes, numPaths] = size(sourcedemand);
jSourcedemand = javaArray('java.lang.Double[]', numTimes);
for iTime = 1:numTimes
    jPath = javaArray('java.lang.Double', numPaths);
    
    for iPath = 1:numPaths
        jPath(iPath) = java.lang.Double(sourcedemand(iTime, iPath));
    end
    jSourcedemand(iTime) = jPath;
end