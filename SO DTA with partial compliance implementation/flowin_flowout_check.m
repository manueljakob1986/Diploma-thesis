function [simulation, T_multiplicator] = flowin_flowout_check(t_mult_start, compliant_demand, noncompliant_demand, ...
    compliant_startvehicles, noncompliant_startvehicles, sink, Yn, Yc)

%[simulation, T_multiplicator] = flowin_flowout_check(t_mult_start, compliant_demand, noncompliant_demand, ...
%    compliant_startvehicles, noncompliant_startvehicles, sink, Yn, Yc)
% flowin_flowout_check compares the total inflow and the total outflow of
% the network and if they are not equal, the time-horizon T is doubled.
%--------------------------------------------------------------------------
% inputs:
% t_mult_start: multiplicator of the time horizon T.
% compliant_demand: input compliant demand matrix (each row is a 
% different time step beginning with t = 0).
% noncompliant_demand: input noncompliant demand matrix (each row is 
% a different time step beginning with t = 0).
% compliant_startvehicles: compliant initial vehicles in each cell i.
% noncompliant_startvehicles: noncompliant initial vehicles in each cell i.
% sink: sink cell.
% Yn: a matrix, showing the number of non-compliant vehicles moving from 
% cell i to cell j at time interval t
% Yc is a matrix, showing the number of compliant vehicles moving from 
% cell i to cell j at time interval t
%
% by Manuel Jakob
% 21 February 2012
%==========================================================================

%total flow-in the network 
%we determine the sum of all initial cars in all cells and the demand
%vehicles for all times
total_flowin = int32(sum(sum(compliant_demand)) + sum(sum(noncompliant_demand)) + ...
    sum(compliant_startvehicles) + sum(noncompliant_startvehicles));

%total flow-out the network
%we determine the sum of all cars in the sink cell for all times
total_flowout = int32(sum(sum(Yn(:,:,sink))) + sum(sum(Yc(:,:,sink))));

T_multiplicator = t_mult_start;
simulation = 0;
if total_flowin ~= total_flowout
    reply = input('The total inflow is not equal to the total outflow of the network. \n The time horizon is set to the double value. \n Do you wish to restart the simulation? \n Yes [Y] \n No  [N] \n (default value [Y]): ', 's');
    if isempty(reply)
        reply = 'Y';
    end
    if reply == 'Y'
        simulation = 1;
        T_multiplicator = T_multiplicator*2;
    elseif reply == 'N'
        simulation = 0;
    end
end
