function [x_veh] = CTM_network_load(dt,T,sourcedemand,paths,discretesplits_nodes)

%Here we are going to load the network based on the
%cell-transmission-model (CTM) from Daganzo. The CTM based network load 
%can be found in the paper from Lo and Szeto. Here we use the TOPL
%simulation tool Sirius - Aurora and the file aurora.jar.
%[x_veh] = CTM_network_load(dt,T,sourcedemand,paths,discretesplits_nodes).
%--------------------------------------------------------------------------
%inputs:
%dt: time interval for the simulation.
%T: total time horizon.
%sourcedemand: demand of the source cell for all paths.
%paths is a (1 x paths)-vector, that shows all existing paths.
%discretesplits_nodes: (2 x number of cells) gives us a correspondance 
%from each link to his beginning and ending node.
%
% by Manuel Jakob
% 11 April 2012
%==========================================================================

%Import class Runner in package runner (Java).
import runner.Runner

%Determine the corresponding turning ratio matrix for all paths
%for all nodes (with an input and output link) at all times.
[paths_TR] = paths_turning_ratios(paths,discretesplits_nodes);
%Translate paths_TR for Java input.
[jSplitRatioMatrix] = java_three_translate(paths_TR);

%Translate sourcedemand for Java input.
[jSourcedemand] = java_two_translate(sourcedemand);

%Run the TOPL simulation tool Sirius - Aurora in Java:
x_veh = Runner.getDensities(dt,dt*T,jSourcedemand,jSplitRatioMatrix);

end











