function [time_dependent_matrix] = time_dependent_changes(fileName, dt)

% [time_dependent_matrix] = time_dependent_changes(fileName, n_time)
% time_dependent_changes reads the file fileName (.xml) and saves all special
% time-dependent values of the fundamental diagram in the matrix 
% time_dependent_matrix. The input dt is the time interval for a 
% number of time-steps due to a given time horizon.
%
% by Manuel Jakob
% 16 February 2012
%==========================================================================

xDoc = xmlread(fullfile(fileName));
xDocRoot = xDoc.getDocumentElement;
scenario = xDoc.getElementsByTagName('scenario').item(0);

%determine the time-dependent changes of the values from the fundamental diagram
profile_set = scenario.getElementsByTagName('fds');
profile_set = profile_set.item(0);
if isempty(profile_set)
    time_dependent_matrix = [];
    return
end
link = profile_set.getElementsByTagName('link');
n_link = link.getLength;
link_id = zeros(1,n_link);
timestep = zeros(1,n_link);
delta = zeros(1,n_link);
qmax = zeros(1,n_link);
N_max = zeros(1,n_link);
time_dependent_matrix = zeros(n_link,5);

for i = 0:n_link-1
        item = link.item(i);
        link_id(i+1) = str2double(item.getAttribute('link_id'));
        timestep(i+1) = str2double(item.getAttribute('timestep'));
        fd = item.getElementsByTagName('fd');
        fd = fd.item(0);
        densityCritical = str2double(fd.getAttribute('densityCritical'));
        densityJam = str2double(fd.getAttribute('densityJam'));
        delta(i+1) = (densityCritical /(densityJam - densityCritical));
        qmax(i+1) = str2double(fd.getAttribute('flowMax'));
        N_max(i+1) = densityJam;

end
        time_dependent_matrix(:,1) = link_id;
        time_dependent_matrix(:,2) = timestep;
        time_dependent_matrix(:,3) = delta;
        time_dependent_matrix(:,4) = qmax/3600 * str2double(dt);
        time_dependent_matrix(:,5) = N_max;
end   
        
       