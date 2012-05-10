function [discretesplits,discretesplits_nodes,add_entry,speed,length] = splitlink(fileName)

% [discretesplits,discretesplits_nodes,add_entry,speed,length] =
% splitlink(fileName)
% splitlink reads the file fileName (.xml) of the network and
% the function output vector discretesplits determines the necessary splits
% of each link. 
% add_entry as the other function output vector gives a relation between 
% all the additional entries (rows or columns) and the splitting links.
% discretesplits_nodes gives us a correspondance from each link to his
% beginning and ending node.
% speed is the free flow speed/speed limit v on the freeways for each link.
% length is the length of each link in the network.
%
% by Manuel Jakob
% 14 February 2012
%==========================================================================

xDoc = xmlread(fullfile(fileName));
xDocRoot = xDoc.getDocumentElement;
scenario = xDoc.getElementsByTagName('scenario').item(0);

%determine the time interval dt
profile_set = scenario.getElementsByTagName('settings');
profile_set = profile_set.item(0);
display = profile_set.getElementsByTagName('display');
display = display.item(0);
dt = str2double(display.getAttribute('dt'));

%determine the free flow speed/speed limit v on the freeways
profile_set = scenario.getElementsByTagName('LinkList');
profile_set = profile_set.item(0);
link = profile_set.getElementsByTagName('link');
n_link = link.getLength;
speed = zeros(1,n_link);
for i = 0:n_link-1
    item = link.item(i);
    fd = item.getElementsByTagName('fd');
    fd = fd.item(0);
    %get the maximal flow for each link
    qmax = str2double(fd.getAttribute('flowMax'));
    %get the critical density for each link
    densityCritical = str2double(fd.getAttribute('densityCritical'));
    %free flow speed/speed limit v on the freeways for each link
    speed(i+1) = qmax/densityCritical;
end    

%(.xml) file of the network part 'LinkList'
profile_set = scenario.getElementsByTagName('LinkList');
profile_set = profile_set.item(0);
link = profile_set.getElementsByTagName('link');
n_link = link.getLength;
length = zeros(1,n_link);
discretesplits = zeros(1,n_link);
discretesplits_nodes = zeros(2,n_link);
add_entry = 0;
count = 0;
for i = 0:n_link-1
    item = link.item(i);
    length(i+1) = str2double(item.getAttribute('length'));
    %Convert the length from miles to feet and 
    %the speed from miles per hour to feets per second.
    discretesplits(i+1) = round(length(i+1)*5280/(dt*speed(i+1)*1.4666666666666666));
    %Look for the beginning node of each link.
    begin = item.getElementsByTagName('begin');
    begin = begin.item(0);
    node_begin = str2double(begin.getAttribute('node_id'));
    %Look for the ending node of each link.
    n_end = item.getElementsByTagName('end');
    n_end = n_end.item(0);
    node_end = str2double(n_end.getAttribute('node_id'));
    %Save beginning and ending node for each link.
    discretesplits_nodes(1,i+1) = - node_begin;
    discretesplits_nodes(2,i+1) = - node_end;
    %add_entry as a function output gives a relation between all the 
    %additional entries (rows or columns) and the splitting links.
    if discretesplits(i+1) > 1
        for j = 1:discretesplits(i+1)-1
            add_entry(count + j) = str2double(item.getAttribute('id'));
        end    
        count = count + j;
    end
end    

