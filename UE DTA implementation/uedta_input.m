% Here we run the function uedta.m with the inputs initial paths and 
% initial flow and the input created by the function 
% [] = topl_myxmlread(fileName, T_multiplicator).
% We also need to specifiy a (.xml) file, stopping criteria and a
% maximum number of iterations.
% The function uedta.m solves the UE DTA problem.
%
% by Manuel Jakob
% 12 March 2012
%==========================================================================

clc; clear all; close all;

%Import the TOPL simulation tool Sirius - Aurora in Java:
sirius_env;

%name of the (.xml) file.
xml_file = 'mySplitSymmetric_complex_aurora.xml';
%xml_file = 'mySplitSymmetric_LongerRightPath_aurora.xml';
%xml_file = 'mySplitSymmetric_aurora.xml';
%xml_file = 'SimpleNetwork.xml';
%define stopping criteria epsilon.
epsilon = 1; %in percent
%define the maximum number of iterations.
iteration_max = 20;
%directory of the (.xml) network file for Java:
xml_file_Java = 'C:\Users\Logge\workspace\MJRunner\data\mySplitSymmetric_complex_sirius.xml';
%xml_file_Java = 'C:\Users\Logge\workspace\MJRunner\data\mySplitSymmetric_LongerRightPath_sirius.xml';
%xml_file_Java = 'C:\Users\Logge\workspace\MJRunner\data\mySplitSymmetric_sirius.xml';
%xml_file_Java = 'C:\Users\Logge\workspace\MJRunner\data\network.xml';

%Run the UE DTA iterative algorithm.
[flow_result,accuracy,accuracy1,accuracy2,flow,iteration,act_tt_result,num_demands] ...
    = uedta(xml_file,epsilon,iteration_max,xml_file_Java);
myXMLEditorDelete(xml_file_Java);

%Plot results for flow on all path vs. number of iterations for the demand existing time steps.
for t = 1:num_demands
    figure
    plot(linspace(0,iteration-1,iteration),flow_result(t,:,1),'r-',linspace(0,iteration-1,iteration),flow_result(t,:,2),'b-',...
        linspace(0,iteration-1,iteration),flow_result(t,:,3),'g-')
    %plot(linspace(0,iteration-1,iteration),flow_result(t,:,1),'r-',linspace(0,iteration-1,iteration),flow_result(t,:,2),'b-')
    xlabel('number of iterations');
    ylabel('flow');
    legend('Path 1','Path 2','Path 3')
    %legend('Shorter Path','Longer Path')
end
 
%Plot accuracy of the convergence criteria vs. number of iterations.
figure
plot(linspace(0,iteration-1,iteration),accuracy,'r-',linspace(0,iteration-1,iteration),accuracy1,'b-',...
    linspace(0,iteration-1,iteration),accuracy2,'g-')
xlabel('number of iterations');
ylabel('accuracy of the convergence criteria');
legend('Accuracy for both time steps','Accuracy time step 1','Accuracy time step 2')

%Plot actual travel times for all paths for the demand existing time steps.
for t = 1:num_demands
    figure
    plot(linspace(0,iteration-1,iteration),act_tt_result(t,:,1),'r-',linspace(0,iteration-1,iteration),...
        act_tt_result(t,:,2),'b-',linspace(0,iteration-1,iteration),act_tt_result(t,:,3),'g-')
    %plot(linspace(0,iteration-1,iteration),act_tt_result(t,:,1),'r-',linspace(0,iteration-1,iteration),act_tt_result(t,:,2))
    xlabel('number of iterations');
    ylabel('actual travel time');
    legend('Path 1','Path 2','Path 3')
    %legend('Shorter Path','Longer Path')
end

