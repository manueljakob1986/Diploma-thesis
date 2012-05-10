function [turning_ratio] = turningratio(demand, node, discretesplits, add_entry, n_link)

%Determine the turning ratios from the (.xml) file property 'knob' in
%'NonCompliantDemandProfileSet'.
%[turning_ratio] = turningratio(demand, node, discretesplits, add_entry, n_link)
%
%--------------------------------------------------------------------------
%Remark: Only a change in the 'knob' value in the first entry of 
%'NonCompliantDemandProfileSet' is valid and it is going to affects all  
%non-compliant flow. There is no fixed turning ratio for compliant flow.
%We get this turning ratio implicitly by solving the Ziliaskopoulos 
%SO DTA LP problem.
%--------------------------------------------------------------------------
%
%inputs:
%demand: demand = profile_set.getElementsByTagName('demand') from the
%(.xml) file part 'NonCompliantDemandProfileSet'.
%node: node = profile_set.getElementsByTagName('node') from the
%(.xml) file part 'NodeList'.
%discretesplits and add_entry are outputs of the function [discretesplits,
%add_entry] = splitlink(fileName).
%n_link: number of links.
%
% by Manuel Jakob
% 11 March 2012
%==========================================================================

%Get the turning ratio from the 'knob' value in the 'NonCompliantDemandProfileSet'
%of the (.xml) file
item = demand.item(0);
tr = str2double(item.getAttribute('knob'));

%Determine the input link of the splitting node.
n_node = node.getLength;
for i = 0:n_node-1 %all nodes
    item = node.item(i);
    outputs = item.getElementsByTagName('outputs');
    outputs = outputs.item(0);
    output = outputs.getElementsByTagName('output');
    n_output = output.getLength;
    
    if n_output > 1
        inputs = item.getElementsByTagName('inputs');
        inputs = inputs.item(0);
        input = inputs.getElementsByTagName('input');
        input = input.item(0);
        input_link = str2double(input.getAttribute('link_id'));
        break
    end
end %all nodes

%Check whether the input link of the splitting node is splitted.
%In case of splitting, change input link to the corresponding splitting
%input link.
if discretesplits(-input_link) > 1
   for i = 1:size(add_entry)
      if add_entry(i) == input_link
          tr_index = i;
      end
   end
else
    tr_index = -input_link;
end

%Calculate the vector turning_ratio
turning_ratio = zeros(1,n_link);
turning_ratio(tr_index) = tr;


