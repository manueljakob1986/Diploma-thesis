function [n_time, num_cells, inc_matrix, q_max_matrix, delta, N_max_matrix, compliant_demand_matrix, ...
    noncompliant_demand_matrix, compliant_initial_vehicles, noncompliant_initial_vehicles, source, sink, dt] ...
    = topl_myxmlread_nosplit(fileName, T_multiplicator)

% [n_time, num_cells, inc_matrix, q_max_matrix, delta, N_max_matrix, compliant_demand_matrix, ...
%    noncompliant_demand_matrix, compliant_initial_vehicles, noncompliant_initial_vehicles, source, sink, dt] = 
%    topl_myxmlread_nosplit(fileName, T_multiplicator)
% topl_myxmlread_nosplit reads the file fileName1 (.xml) of the network
% (no splitting of the links) and has the following outputs:
%--------------------------------------------------------------------------
% n_time: total time horizon.
% num_cells: number of cells after splitting some links that are too long.
% inc_matrix: incidence matrix of the network.
% q_max_matrix: maximum number of vehicles that can flow in or out cell i 
% during the t-th time period (each row is a different time step).
% delta: free-flow to backward propagation speed ratio for cell i at time t 
% (for each time step t the diagonal of the matrix delta(t) changes).
% N_max_matrix: maximum number of vehicles in cell i at time t.
% compliant_demand_matrix: input compliant demand matrix (each row is a 
% different time step beginning with t = 0).
% noncompliant_demand_matrix: input noncompliant demand matrix (each row is 
% a different time step beginning with t = 0).
% turning_ratio: vector of dimension num_cells including all turning ratios
% for all splitting cells.
% compliant_initial_vehicles: compliant initial vehicles in each cell i.
% noncompliant_initial_vehicles: noncompliant initial vehicles in each cell i.
% source: source cell.
% sink: sink cell.
% dt: time interval for a number of time-steps due to a given time horizon.
%
% by Manuel Jakob
% 06 February 2012
%==========================================================================

xDoc = xmlread(fullfile(fileName));
xDocRoot = xDoc.getDocumentElement;
scenario = xDoc.getElementsByTagName('scenario').item(0);

%==========================================================================
%determine the time horizon n_time
profile_set = scenario.getElementsByTagName('settings');
profile_set = profile_set.item(0);
display = profile_set.getElementsByTagName('display');
display = display.item(0);
dt = str2double(display.getAttribute('dt'));
timeMax = display.getAttribute('timeMax');
n_time = str2double(timeMax)/dt;
n_time = n_time*T_multiplicator;

%==========================================================================
%(.xml) file part 'LinkList'
profile_set = scenario.getElementsByTagName('LinkList');
profile_set = profile_set.item(0);
link = profile_set.getElementsByTagName('link');
n_link = link.getLength;
    
%initialize variables    
q_max_matrix = zeros(n_time,n_link);
N_max_matrix = zeros(n_time,n_link);
delta = zeros(n_time,n_link);
length = zeros(1,n_link);
for i = 0:n_link-1
    item = link.item(i);

    %create a vector of the lengths of all other links from the .xml file
    length(i+1) = str2double(item.getAttribute('length'));
             
    %create the matrix for q_max of all other links from the .xml file
    fd = item.getElementsByTagName('fd');
    fd = fd.item(0);
    qmax = fd.getAttribute('flowMax');
    q_max_matrix(:,i+1) = str2double(qmax)/3600 * dt * ones(n_time,1);
          
    %create the free-flow to backward propagation speed ratio delta
    %of all other links from the .xml file
    densityCritical = fd.getAttribute('densityCritical');
    densityJam = fd.getAttribute('densityJam');
    delta(:,i+1) = (str2double(densityCritical) /(str2double(densityJam) - ...
        str2double(densityCritical))) * ones(n_time,1);
        
    %create the matrix for N_max of all other links from the .xml file
    N_max = densityJam;
    N_max_matrix(:,i+1) = str2double(N_max) * ones(n_time,1);          
end
num_cells = n_link;


%==========================================================================
% Change all special time-dependent values of the fundamental diagram due to the file 
% extras.xml and the function [time_dependent_matrix] = time_dependent_changes(fileName).
% All other values for each link for different time-steps that are not 
% mentionned in file extras.xml are set as the standard values for each
% link from the network (.xml) file.

[time_dependent_matrix] = time_dependent_changes(fileName, dt);

for i = 1:size(time_dependent_matrix,1)             
    %change the time-dependent value for all the following
    %time-steps in the matrices q_max_matrix and delta
    for k = time_dependent_matrix(i,2):n_time %change for all following time-steps 
        q_max_matrix(k,-time_dependent_matrix(i,1)) = time_dependent_matrix(i,4);
        delta(k,-time_dependent_matrix(i,1)) = time_dependent_matrix(i,3);
        N_max_matrix(k,-time_dependent_matrix(i,1)) = time_dependent_matrix(i,5);
    end %change for all following time-steps 
end

   
%==========================================================================
%(.xml) file part 'CompliantInitialDensityProfile'
profile_set = scenario.getElementsByTagName('CompliantInitialDensityProfile');
profile_set = profile_set.item(0);
density = profile_set.getElementsByTagName('density');
n_compliant_density = density.getLength;
%create the vector compliant_initial_vehicles from the .xml file
compliant_initial_vehicles = zeros(1,num_cells);
for i = 0:n_compliant_density-1
    item = density.item(i);
    
    %create a vector due to the compliant intial vehicles for the splitting links from the .xml file
    compliant_initial_vehicles(i+1) = str2double(item.getTextContent) * length(i+1);
end

    
%==========================================================================
%(.xml) file part 'NonCompliantInitialDensityProfile'
profile_set = scenario.getElementsByTagName('NonCompliantInitialDensityProfile');
profile_set = profile_set.item(0);
density = profile_set.getElementsByTagName('density');
n_noncompliant_density = density.getLength;
%create the vector noncompliant_initial_vehicles from the .xml file
noncompliant_initial_vehicles = zeros(1,num_cells);

for i = 0:n_noncompliant_density-1
    item = density.item(i);
    %create a vector due to the noncompliant intial vehicles for the splitting links from the .xml file
    noncompliant_initial_vehicles(i+1) = str2double(item.getTextContent) * length(i+1);
end
    

%==========================================================================    
%create the incidence matrix inc_matrix from the .xml file
inc_matrix = zeros(n_link,n_link);
profile_set = scenario.getElementsByTagName('NodeList');
profile_set = profile_set.item(0);
node = profile_set.getElementsByTagName('node');
n_node = node.getLength;
for i = 0:n_node-1 %all nodes
    item = node.item(i);
    %look for input links
    inputs = item.getElementsByTagName('inputs');
    inputs = inputs.item(0);
    input = inputs.getElementsByTagName('input');
    n_input = input.getLength;
    %determine the source cell
    if n_input == 0
        outputs = item.getElementsByTagName('outputs');
        outputs = outputs.item(0);
        output = outputs.getElementsByTagName('output');
        output = output.item(0);
        source = - str2double(output.getAttribute('link_id'));
    end
    
    for j = 0:n_input-1 %all input links
        %save all input links for each node
        input = inputs.getElementsByTagName('input');
        input = input.item(j); 
        %input link
        a = str2double(input.getAttribute('link_id')); 
        
        %look for output links
        outputs = item.getElementsByTagName('outputs');
        outputs = outputs.item(0);
        output = outputs.getElementsByTagName('output');
        n_output = output.getLength; 
        %determine the sink cell
        if n_output == 0
            inputs = item.getElementsByTagName('inputs');
            inputs = inputs.item(0);
            input = inputs.getElementsByTagName('input');
            input = input.item(0);
            sink = - str2double(input.getAttribute('link_id'));
        end 
        for k = 0:n_output-1 %all output links
            output = outputs.getElementsByTagName('output');
            output = output.item(k);
            b = str2double(output.getAttribute('link_id'));   
            %save all output links for each node and the corresponding input
            inc_matrix(-a,-b) = 1;
        end %all output links
   end %all intput links
end %all nodes 


%==========================================================================    
%(.xml) file part 'CompliantDemandProfileSet'
profile_set = scenario.getElementsByTagName('CompliantDemandProfileSet');
profile_set = profile_set.item(0);
demand = profile_set.getElementsByTagName('demand');
n_compliant_demand_links = demand.getLength;
compliant_demand_matrix = zeros(n_time,num_cells);

for j = 0:n_compliant_demand_links-1 %all demand links
    item = demand.item(j);
    linkid = str2double(item.getAttribute('link_id'));
    splitd = strsplit(char(item.getTextContent),', ');
    n_compliant_demand = size(splitd,2);

    %create the matrix compliant_demand_matrix from the .xml file,
    %different time-steps for all not splitted demand cells
    for i = 0:n_compliant_demand-1
        compliant_demand_matrix(i+1,-linkid) = str2double(splitd{i+1});
        %last time-dependent demand value is set equal to
        %all following demand values in compliant_demand_matrix
        if i == n_compliant_demand-1
            for k = i+2:size(compliant_demand_matrix,1)
                compliant_demand_matrix(k,-linkid) = str2double(splitd{i+1});
            end
        end
    end
end %all demand links


%==========================================================================    
%(.xml) file part 'NonCompliantDemandProfileSet'
profile_set = scenario.getElementsByTagName('NonCompliantDemandProfileSet');
profile_set = profile_set.item(0);
demand = profile_set.getElementsByTagName('demand');
n_noncompliant_demand_links = demand.getLength;
noncompliant_demand_matrix = zeros(n_time,num_cells);

for j = 0:n_noncompliant_demand_links-1 %all demand links
    item = demand.item(j);
    linkid = str2double(item.getAttribute('link_id'));
    splitd = strsplit(char(item.getTextContent),', ');
    n_noncompliant_demand = size(splitd,2);

    %create the matrix noncompliant_demand_matrix from the .xml file,
    %different time-steps for all not splitted demand cells
    for i = 0:n_noncompliant_demand-1
        noncompliant_demand_matrix(i+1,-linkid) = str2double(splitd{i+1});
        %last time-dependent demand value is set equal to
        %all following demand values in noncompliant_demand_matrix
        if i == n_noncompliant_demand-1
            for k = i+2:size(noncompliant_demand_matrix,1)
                noncompliant_demand_matrix(k,-linkid) = str2double(splitd{i+1});
            end
        end
    end 
end %all demand links
