function [n_time, num_cells, inc_matrix, q_max_matrix, delta, N_max_matrix, compliant_demand_matrix, ...
    noncompliant_demand_matrix, turning_ratio, compliant_initial_vehicles, noncompliant_initial_vehicles, source, sink, dt] ...
    = topl_myxmlread(fileName, T_multiplicator)

% [n_time, num_cells, inc_matrix, q_max_matrix, delta, N_max_matrix, compliant_demand_matrix, ...
%    noncompliant_demand_matrix, turning_ratio, compliant_initial_vehicles, noncompliant_initial_vehicles, source, sink, dt] = 
%    topl_myxmlread(fileName, T_multiplicator)
% topl_myxmlread reads the file fileName1 (.xml) of the network and has the
% following outputs:
%--------------------------------------------------------------------------
% n_time: time-steps the a given time-horizon.
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
n_link_old = n_link;

%determine links that we need to split because they are too long
    length = zeros(1,n_link);
    [discretesplits,~,add_entry,~,~] = splitlink(fileName);
    link_id = 0;
    index = 1;
    for i = 0:n_link-1
        item = link.item(i);
        length(i+1) = str2double(item.getAttribute('length'));
        if discretesplits(i+1) > 1
            new_cells = discretesplits(i+1) - 1;
            n_link = n_link + new_cells;
            link_id(index) = str2double(item.getAttribute('id')); 
            index = index + 1;
        end
    end 
    
%initialize variables    
q_max_matrix = zeros(n_time,n_link);
N_max_matrix = zeros(n_time,n_link);
delta = zeros(n_time,n_link);
length = zeros(1,n_link);
columns = 0;
    for i = 0:n_link_old-1
        item = link.item(i);
        
        for index = 1:size(link_id,2)
            % if the link needs to be splitted
            if str2double(item.getAttribute('id')) == link_id(index)
                
                %create a vector of the lengths for the splitting links from the .xml file
                length(i+1) = str2double(item.getAttribute('length'))/discretesplits(i+1);
                %split length of link into discretesplits(i+1) number of parts
                for j = 1:discretesplits(i+1)-1
                    length(n_link_old + j + columns) = str2double(item.getAttribute('length'))/discretesplits(i+1);
                end
                
                %create the matrix for q_max for the splitting links from the .xml file
                fd = item.getElementsByTagName('fd');
                fd = fd.item(0);
                qmax = fd.getAttribute('flowMax');
                q_max_matrix(:,i+1) = str2double(qmax)/3600 * dt * ones(n_time,1);
                %create discretesplits(i+1)-1 new columns in the q_max_matrix
                for j = 1:discretesplits(i+1)-1
                    q_max_matrix(:,n_link_old + j + columns) = str2double(qmax)/3600 * dt * ones(n_time,1);
                end
                    
                %create the free-flow to backward propagation speed ratio delta
                %for the splitting links from the .xml file
                densityCritical = fd.getAttribute('densityCritical');
                densityJam = fd.getAttribute('densityJam');
                delta(:,i+1) = (str2double(densityCritical) /(str2double(densityJam) - ...
                    str2double(densityCritical))) * ones(n_time,1);
                %create discretesplits(i+1)-1 new columns in the delta matrix
                for j = 1:discretesplits(i+1)-1
                delta(:,n_link_old + j + columns) = (str2double(densityCritical) /(str2double(densityJam) - ...
                    str2double(densityCritical))) * ones(n_time,1);
                end
        
                %create the matrix for N_max for the splitting links from the .xml file
                N_max = densityJam;
                N_max_matrix(:,i+1) = str2double(N_max) * ones(n_time,1);
                %create discretesplits(i+1)-1 new columns in the N_max_matrix
                for j = 1:discretesplits(i+1)-1
                    N_max_matrix(:,n_link_old + j + columns) = str2double(N_max) * ones(n_time,1);
                end
                
                %count all extra columns for the splitted links
                columns = columns + discretesplits(i+1) - 1; 
                
            break    
            else %when there is no split of the link
                if index == size(link_id,2)

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
            end
        end
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

    for index = 1:size(link_id,2)
            if time_dependent_matrix(i,1) == link_id(index)
            % if the link needs to be splitted
            
                %change the time-dependent value for all the following
                %time-steps in the matrices q_max_matrix and delta
                for k = time_dependent_matrix(i,2):n_time %change for all following time-steps   
                    q_max_matrix(k,-time_dependent_matrix(i,1)) = time_dependent_matrix(i,4);
                    delta(k,-time_dependent_matrix(i,1)) = time_dependent_matrix(i,3);
                    N_max_matrix(k,-time_dependent_matrix(i,1)) = time_dependent_matrix(i,5);
                
                    for l = 1:size(add_entry,2)
                        if add_entry(l) == time_dependent_matrix(i,1)
                            q_max_matrix(k,n_link_old + l) = time_dependent_matrix(i,4);
                            delta(k,n_link_old + l) = time_dependent_matrix(i,3);
                            N_max_matrix(k,n_link_old + l) = time_dependent_matrix(i,5);
                        end    
                    end
                end %change for all following time-steps
                
            break    
            else
            % if the link does not need to be splitted
                if index == size(link_id,2)
                  
                    %change the time-dependent value for all the following
                    %time-steps in the matrices q_max_matrix and delta
                    for k = time_dependent_matrix(i,2):n_time %change for all following time-steps 
                    q_max_matrix(k,-time_dependent_matrix(i,1)) = time_dependent_matrix(i,4);
                    delta(k,-time_dependent_matrix(i,1)) = time_dependent_matrix(i,3);
                    N_max_matrix(k,-time_dependent_matrix(i,1)) = time_dependent_matrix(i,5);
                    end %change for all following time-steps 
                end
            end
    end       
    
end

   
%==========================================================================
%(.xml) file part 'CompliantInitialDensityProfile'
profile_set = scenario.getElementsByTagName('CompliantInitialDensityProfile');
profile_set = profile_set.item(0);
density = profile_set.getElementsByTagName('density');
n_compliant_density = density.getLength;
%create the vector compliant_initial_vehicles from the .xml file
compliant_initial_vehicles = zeros(1,num_cells);
entries = 0;
    for i = 0:n_compliant_density-1
        item = density.item(i);
        
        for index = 1:size(link_id,2)
            % if the link needs to be splitted
            if str2double(item.getAttribute('link_id')) == link_id(index)
               %initial density is splitted to conserve the original density at the end
               split_density = str2double(item.getTextContent)/discretesplits(i+1);
               %create a vector due to the intial vehicles for the splitting links from the .xml file
               compliant_initial_vehicles(i+1) = split_density * length(i+1);
                %create discretesplits(i+1)-1 new entries in the compliant_initial_vehicles vector
                for j = 1:discretesplits(i+1)-1
                    compliant_initial_vehicles(n_compliant_density + j + entries) = split_density * length(n_compliant_density + j + entries); 
                end
                
                %count all extra entries for the splitted links
                entries = entries + discretesplits(i+1) - 1;
                
            break    
            else %when there is no split of the link
               if index == size(link_id,2) 
                   %create a vector due to the compliant intial vehicles for the splitting links from the .xml file
                   compliant_initial_vehicles(i+1) = str2double(item.getTextContent) * length(i+1);
               end     
            end
        end    
    end

    
%==========================================================================
%(.xml) file part 'NonCompliantInitialDensityProfile'
profile_set = scenario.getElementsByTagName('NonCompliantInitialDensityProfile');
profile_set = profile_set.item(0);
density = profile_set.getElementsByTagName('density');
n_noncompliant_density = density.getLength;
%create the vector noncompliant_initial_vehicles from the .xml file
noncompliant_initial_vehicles = zeros(1,num_cells);
entries = 0;
    for i = 0:n_noncompliant_density-1
        item = density.item(i);
        
        for index = 1:size(link_id,2)
            % if the link needs to be splitted
            if str2double(item.getAttribute('link_id')) == link_id(index)
               %initial density is splitted to conserve the original density at the end
               split_density = str2double(item.getTextContent)/discretesplits(i+1);
               %create a vector due to the intial vehicles for the splitting links from the .xml file
               noncompliant_initial_vehicles(i+1) = split_density * length(i+1);
                %create discretesplits(i+1)-1 new entries in the noncompliant_initial_vehicles vector
                for j = 1:discretesplits(i+1)-1
                    noncompliant_initial_vehicles(n_noncompliant_density + j + entries) = split_density * length(n_noncompliant_density + j + entries); 
                end
                
                %count all extra entries for the splitted links
                entries = entries + discretesplits(i+1) - 1;
                
            break    
            else %when there is no split of the link
               if index == size(link_id,2) 
                   %create a vector due to the noncompliant intial vehicles for the splitting links from the .xml file
                   noncompliant_initial_vehicles(i+1) = str2double(item.getTextContent) * length(i+1);
               end     
            end
        end    
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
        
            for input_index = 1:size(link_id,2)
                if str2double(input.getAttribute('link_id')) == link_id(input_index)
                %----------------------------------------------------------
                %case 1 and 2: if input links needs to splitted,                        
                %for output link splitted and output link non-splitted

                    for k = 0:n_output-1 %all output links
                        output = outputs.getElementsByTagName('output');
                        output = output.item(k);
                        for output_index = 1:size(link_id,2)
                            % if the output link needs to be splitted
                            if str2double(output.getAttribute('link_id')) == link_id(output_index)
                            b = str2double(output.getAttribute('link_id')); 
                        
                            %count the number of splitting partial links
                            %for a link a that needs to be splitted
                            count_id_input = 0;
                            for l = 1:size(add_entry,2)
                                if add_entry(l) == a
                                    count_id_input = count_id_input + 1;
                                end    
                            end 
                            %look for the links b in the vector add_entry
                            %and save all output links for each node and the
                            %corresponding input
                            for l = 1:size(add_entry,2)
                                if add_entry(l) == b
                                    inc_matrix(-a,n_link_old + l) = 1;
                                        for m = (l+1):size(add_entry,2)
                                            if add_entry(m) == b 
                                                inc_matrix(n_link_old + m - 1,n_link_old + m) = 1;
                                            else
                                                l = m - 1;
                                            break
                                            end
                                        end
                                    inc_matrix(n_link_old + l,-b) = 1;
                                    break
                                end    
                            end
                            %look for the link a in the vector add_entry
                            %and save all output links for each node and the
                            %corresponding input
                            count_link = 1;
                            for l = 1:size(add_entry,2)
                                if add_entry(l) == a
                                        for m = (l+1):size(add_entry,2)
                                            if count_link < count_id_input
                                                inc_matrix(n_link_old + m - 1,n_link_old + m) = 1;
                                            elseif count_link == count_id_input
                                                l = m - 1;
                                            break
                                            end
                                            count_link = count_link + 1;
                                        end
                                    inc_matrix(n_link_old + l,-a) = 1;
                                    break
                                end    
                            end
                            break
                            
                            %==============================================
                            else %output links do not need to be splitted
                                if output_index == size(link_id,2)
                                    b = str2double(output.getAttribute('link_id')); 
                        
                                    %count the number of splitting partial links
                                    %for link a that needs to be splitted
                                    count_id = 0;
                                    for l = 1:size(add_entry,2)
                                        if add_entry(l) == a
                                            count_id = count_id + 1;
                                        end    
                                    end 
                                    %look for the link a in the vector add_entry
                                    %and save all output links for each node and the
                                    %corresponding input
                                    count_link = 1;
                                    for l = 1:size(add_entry,2)
                                        if add_entry(l) == a
                                            inc_matrix(-a,-b) = 1;
                                                for m = (l+1):size(add_entry,2)
                                                    if count_link < count_id
                                                        inc_matrix(n_link_old + m - 1,n_link_old + m) = 1;
                                                    elseif count_link == count_id
                                                        l = m - 1;
                                                    break
                                                    end
                                                    count_link = count_link + 1;
                                                end
                                            inc_matrix(n_link_old + l,-a) = 1;
                                            break
                                        end    
                                    end 
                                end % if output_index == size(link_id,2)
                            end % if the output link needs to be splitted
                        end %for output_index                     
                    end %all output links
                
                break    
                else
                if input_index == size(link_id,2)    
                %----------------------------------------------------------
                %case 3 and 4: if input links do not need to splitted,                        
                %for output link splitted and output link non-splitted
        
                    for k = 0:n_output-1 %all output links
                        output = outputs.getElementsByTagName('output');
                        output = output.item(k);
                        for output_index = 1:size(link_id,2)
                            % if the output link needs to be splitted
                            if str2double(output.getAttribute('link_id')) == link_id(output_index)
                            b = str2double(output.getAttribute('link_id')); 
                            
                            %look for the link b in the vector add_entry
                            %and save all output links for each node and the
                            %corresponding input
                            for l = 1:size(add_entry,2)
                                if add_entry(l) == b
                                    inc_matrix(-a,n_link_old + l) = 1;
                                        for m = (l+1):size(add_entry,2)
                                            if add_entry(m) == b 
                                                inc_matrix(n_link_old + m - 1,n_link_old + m) = 1;
                                            else
                                                l = m - 1;
                                            break
                                            end
                                        end
                                    inc_matrix(n_link_old + l,-b) = 1;
                                    break
                                end    
                            end
                            break
                            
                            %==============================================
                            else %output links do not need to be splitted
                                if output_index == size(link_id,2)
                                    b = str2double(output.getAttribute('link_id')); 
                        
                                    %save all output links for each node and the
                                    %corresponding input
                                    inc_matrix(-a,-b) = 1;
                                end
                            end % if the output link needs to be splitted
                        end %for output_index
                    end %all output links 
                end %if input_index == size(link_id,2)  
                end %if input links needs to splitted
            end %for input_index     
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
    
    for index = 1:size(link_id,2)
        % if the link needs to be splitted
        if linkid == link_id(index)
            %create the matrix compliant_demand_matrix from the .xml file,
            %different time-steps for the splitted demand cells,
            %demand is splitted to conserve the original demand for
            %splitted links
            for i = 0:n_compliant_demand-1     
                compliant_demand_matrix(i+1,-linkid) = str2double(splitd{i+1})/discretesplits(-linkid);
                    %last time-dependent demand value is set equal to
                    %all following demand values in compliant_demand_matrix
                    if i == n_compliant_demand-1
                        for k = i+2:size(compliant_demand_matrix,1)
                            compliant_demand_matrix(k,-linkid) = str2double(splitd{i+1})/discretesplits(-linkid);
                        end
                    end
            end  
            
            for l = 1:size(add_entry,2)
                if add_entry(l) == linkid
                    for i = 0:n_compliant_demand-1
                        compliant_demand_matrix(i+1,n_link_old + l) = str2double(splitd{i+1})/discretesplits(-linkid);
                            %last time-dependent demand value is set equal to
                            %all following demand values in compliant_demand_matrix
                            if i == n_compliant_demand-1
                                for k = i+2:size(compliant_demand_matrix,1)
                                    compliant_demand_matrix(k,n_link_old + l) = str2double(splitd{i+1})/discretesplits(-linkid);
                                end
                            end
                    end
                end
            end
            
        break    
        else %when there is no split of the link
            if index == size(link_id,2)

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
            end
        end %if splitted
    end %for index = 1:size(link_id,2)
end %all demand links


%==========================================================================    
%(.xml) file part 'NonCompliantDemandProfileSet'
profile_set = scenario.getElementsByTagName('NonCompliantDemandProfileSet');
profile_set = profile_set.item(0);
demand = profile_set.getElementsByTagName('demand');
n_noncompliant_demand_links = demand.getLength;

%Determine the vector turning_ratio with the function: 
%[turning_ratio] = turningratio(demand, node, discretesplits, add_entry, n_link)
[turning_ratio] = turningratio(demand, node, discretesplits, add_entry, n_link);

noncompliant_demand_matrix = zeros(n_time,num_cells);

for j = 0:n_noncompliant_demand_links-1 %all demand links
    item = demand.item(j);
    linkid = str2double(item.getAttribute('link_id'));
    splitd = strsplit(char(item.getTextContent),', ');
    n_noncompliant_demand = size(splitd,2);
    
    for index = 1:size(link_id,2)
        % if the link needs to be splitted
        if linkid == link_id(index)
            %create the matrix noncompliant_demand_matrix from the .xml file,
            %different time-steps for the splitted demand cells,
            %demand is splitted to conserve the original demand for
            %splitted links
            for i = 0:n_noncompliant_demand-1     
                noncompliant_demand_matrix(i+1,-linkid) = str2double(splitd{i+1})/discretesplits(-linkid);
                    %last time-dependent demand value is set equal to
                    %all following demand values in noncompliant_demand_matrix
                    if i == n_noncompliant_demand-1
                        for k = i+2:size(noncompliant_demand_matrix,1)
                            noncompliant_demand_matrix(k,-linkid) = str2double(splitd{i+1})/discretesplits(-linkid);
                        end
                    end
            end  
            
            for l = 1:size(add_entry,2)
                if add_entry(l) == linkid
                    for i = 0:n_noncompliant_demand-1
                        noncompliant_demand_matrix(i+1,n_link_old + l) = str2double(splitd{i+1})/discretesplits(-linkid);
                            %last time-dependent demand value is set equal to
                            %all following demand values in noncompliant_demand_matrix
                            if i == n_noncompliant_demand-1
                                for k = i+2:size(noncompliant_demand_matrix,1)
                                    noncompliant_demand_matrix(k,n_link_old + l) = str2double(splitd{i+1})/discretesplits(-linkid);
                                end
                            end
                    end
                end
            end
            
        break    
        else %when there is no split of the link
            if index == size(link_id,2)

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
            end
        end %if splitted
    end %for index = 1:size(link_id,2)
end %all demand links
