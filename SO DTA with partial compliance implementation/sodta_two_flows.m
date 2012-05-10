function [simulation, T_multiplicator] = sodta_two_flows(T, numcells, sink, incidencematrix, maxvehflow, speedratio, ...
    maxvehin, compliant_demand, noncompliant_demand, turning_ratios, compliant_startvehicles, noncompliant_startvehicles, t_mult_start)
%The function sodta_two_flows is a modified version of the Matlab function
%sodta that is extended to the Ziliaskopoulos framework to handle two types
%of flows. These are compliant and non-compliant flows. The non-compliant
%flow will have fixed turning ratios and not be controlled by the
%optimization problem. The compliant flow is the flow that we can control and
%the flow whose travel time we wish to minimize.
%--------------------------------------------------------------------------
%inputs:
%demand is the input demand matrix (each row is a different time step beginning with t = 0).
%maxvehflow is the maximum number of vehicles that can flow in or out cell i during the t-th time period 
%(each row is a different time step).
%delta(i)=diag(speedratio(i,:)) is the free-flow to backward propagation speed ratio for cell i at time t 
%(for each time step t the diagonal of the matrix delta(t) changes).
%maxvehin is the maximum number of vehicles in cell i at time t.
%startvehicles is a vector with dimension numcells.
%turning_ratios is a vector of dimension numcells including all turning ratios for all cells. Every entry of 
%turning_ratios is zero when there is no split, otherwise the entry is the splitting ratio of the junction.
%
% by Manuel Jakob
% 29 November 2011
%==========================================================================


%cases that are not possible:
if (size(maxvehflow,1)~=T)
    display('This is not a valid network')
    return
end
if (size(speedratio,1)~=T)
    display('This is not a valid network')
    return
end
if (size(maxvehin,1)~=T)
    display('This is not a valid network')
    return
end
if (size(compliant_demand,1)~=T)
    display('This is not a valid network')
    return
end
if (size(noncompliant_demand,1)~=T)
    display('This is not a valid network')
    return
end
if (size(maxvehflow,2)~=numcells)
    display('This is not a valid network')
    return
end
if (size(speedratio,2)~=numcells)
    display('This is not a valid network')
    return
end
if (size(maxvehin,2)~=numcells)
    display('This is not a valid network')
    return
end
if (size(compliant_demand,2)~=numcells)
    display('This is not a valid network')
    return
end
if (size(noncompliant_demand,2)~=numcells)
    display('This is not a valid network')
    return
end
if (size(turning_ratios,2)~=numcells)
    display('This is not a valid network')
    return
end
if (size(compliant_startvehicles,2)~=numcells)
    display('This is not a valid network')
    return
end
if (size(noncompliant_startvehicles,2)~=numcells)
    display('This is not a valid network')
    return
end

%==========================================================================
cvx_begin

    variable xc(T+1,numcells);                %xc is a matrix, showing the number of compliant vehicles 
                                              %in cell i at time interval t
    variable xn(T+1,numcells);                %xn is a matrix, showing the number of non-compliant vehicles 
                                              %in cell i at time interval t
    variable Yc(T+1,numcells,numcells);       %Yc is a matrix, showing the number of compliant vehicles 
                                              %moving from cell i to cell j at time interval t
    variable Yn(T+1,numcells,numcells);       %Yn is a matrix, showing the number of non-compliant vehicles 
                                              %moving from cell i to cell j at time interval t

%define the matrix Yc and Yn (with the help of the incidence matrix) such that Yc(t) = Yn(t) = 0, 
%whenever two vertices are not connected
for i=1:numcells
    for j=1:numcells
        if incidencematrix(i,j) == 0
            Yc(:,i,j) == 0;
            Yn(:,i,j) == 0; 
        end
    end
end

%define vector e 
e = zeros(1, numcells);
for i=1:numcells
    e(i) = 1;
    if i == sink
        e(i) = 0;
    end
end

%define objective vector c
c = [];
for i=1:T
    c =  [c,e];
end

%initial conditions
Yc(1,:,:) == zeros(1,numcells,numcells);
Yn(1,:,:) == zeros(1,numcells,numcells);

xc(1,:) == compliant_startvehicles;
xn(1,:) == noncompliant_startvehicles;
for j=1:T
    xc(j+1,:) == xc(j,:) + (((reshape(Yc(j,:,:),numcells,numcells))' - ...
        reshape(Yc(j,:,:),numcells,numcells))*ones(numcells,1))' + compliant_demand(j,:);
    xn(j+1,:) == xn(j,:) + (((reshape(Yn(j,:,:),numcells,numcells))' - ...
        reshape(Yn(j,:,:),numcells,numcells))*ones(numcells,1))' + noncompliant_demand(j,:);
end

%define optimization vector zc for compliant vehicles and zn for
%non-compliant vehicles
zc = [];
zn = [];
for i=2:T+1
    zc = [zc;(xc(i,:))'];
    zn = [zn;(xn(i,:))'];
end

%define matrix Ac for compliant vehicles and An for non-compliant for all times T 
%for restrictions with Yc(t) and Yn(t)
Ac = [];
An = [];
for i=2:T+1
    Ac = [Ac ; reshape(Yc(i,:,:),numcells,numcells)];
    An = [An ; reshape(Yn(i,:,:),numcells,numcells)];
end
%define matrix Bc for compliant vehicles and Bn for non-compliant for all times T 
%for restrictions with Yc(t)^T and Yn(t)^T
Bc = [];
Bn = [];
for i=2:T+1
    Bc = [Bc ; (reshape(Yc(i,:,:),numcells,numcells))'];
    Bn = [Bn ; (reshape(Yn(i,:,:),numcells,numcells))'];
end
%define the vector Q for maxvehflow for all times T
Q = [];
for i=1:T
    Q = [Q ; (maxvehflow(i,:))'];
end
%define the vector N for maxvehin for all times T
N = [];
for i=1:T
    N = [N ; (maxvehin(i,:))'];
end
%define the block-diagonal matrix delta for the speed ratio for all times T
delta = [];
for i=1:T
    delta = blkdiag(delta, diag(speedratio(i,:)));
end
%define the diagonal matrix TR for the turning ratios for all times T
TR = [];
for i=1:T
    TR = blkdiag(TR, diag(turning_ratios));
end

%define a help vector h for the split ratio constraint for non-compliant
%users with the help of the vector turning_ratios
h = zeros(length(turning_ratios),1);
for i = 1:length(turning_ratios)
     if turning_ratios(i) ~= 0
         for j = 1:numcells
             if incidencematrix(i,j) ~= 0
             h(j) = 1;
             break
             end
         end    
     end
end

%Creating a help vector for the additional constraint formulation to reach 
%User Equilibrium state.
%svector = sink-1:numcells:(sink-1 + (T-1)*numcells);
%Creating the penalty part of the objective function:
penalty = 0;
for t = 2:T+1
    for i = 1:numcells
        for j = 1:numcells
            penalty = penalty + (Yc(t,i,j) + Yn(t,i,j))*(T + 1 - t);
        end
    end
end
%Solve the linear Programm
%minimize total travel time
minimize(c*zc + c*zn);
%minimize total travel time with penalty function to suppress holding
%minimize((c*zc + c*zn)*(T+1) - penalty);
subject to 
    zc >= 0;
    zn >= 0;
    Ac >= 0;
    An >= 0;
    Ac*ones(numcells,1) + An*ones(numcells,1) <= zc + zn;
    Ac*ones(numcells,1) + An*ones(numcells,1) <= Q;
    Bc*ones(numcells,1) + Bn*ones(numcells,1) <= Q;
    Bc*ones(numcells,1) + Bn*ones(numcells,1) <= delta*(N - zc - zn);
    An*h == TR*An*ones(numcells,1);
    
    Ac*ones(numcells,1) <= zc;
    An*ones(numcells,1) <= zn;
    
    %additional constraint to reach User Equilibrium state:
    for t = 1:T
        Yc(t+1,5,6) + Yn(t+1,5,6) <= speedratio(t,5)*(maxvehin(t,5) - xc(t+1,5) - xn(t+1,5));
    end
%     Yn(2,5,6) <= 35.5556;
%     Yn(3,5,6) <= 35.5556;
%     Yn(4,5,6) <= 35.5556;
%     Yn(5,5,6) <= 35.5556;
%     Yn(6,5,6) <= 35.5556;
%     Ac(svector,:)*ones(numcells,1) + An(svector,:)*ones(numcells,1) <= ...
%         delta(svector,svector)*(N(svector,:) - zc(svector,:) - zn(svector,:));
%     Ac(svector+6,:)*ones(numcells,1) + An(svector+6,:)*ones(numcells,1) <= ...
%         delta(svector+6,svector+6)*(N(svector+6,:) - zc(svector+6,:) - zn(svector+6,:));

cvx_end
%==========================================================================
% for t = 1:T
% test0(t)=    Yn(t+1,4,5);
% test1(t)=    Yn(t+1,5,6);
% test2(t)=    xn(t+1,4);
% test3(t)=    xn(t+1,5);
% test4(t)=speedratio(t,5)*(maxvehin(t,5) - xc(t+1,5) - xn(t+1,5));
% end
% test0
% test1
% test2
% test3
% test4
% xn(1:6,1:6)
% Yn(1:6,1:6,1:6)

%testing the FIFO property (approximate)
%left side of the FIFO property:
LH = (Ac*ones(numcells,1))./(Ac*ones(numcells,1) + An*ones(numcells,1));
%right side of the FIFO property:
RH = zc./(zc + zn);

%determine the approximation error for the given network
error = LH - RH;

%replaces zeros in the denominator of the LH through zeros in the 
%component of error
for i = 1:T*numcells
    if Ac(i,:)*ones(numcells,1) + An(i,:)*ones(numcells,1) == 0
        error(i) = 0;
    end
end

%==========================================================================
%count all nonzero entries in the incidence matrix
count = 0;
for i = 1:numcells
    for j = 1:numcells
        if incidencematrix(i,j) ~= 0
            count = count + 1;
        end
    end
end
        
%plot Yn, Yc only for the nonzero entries for all times T
k = 0;
t = 1:T;
for i = 1:numcells
    for j = 1:numcells
        if incidencematrix(i,j) ~= 0
            k = k + 1;
            subplot(ceil(count/3),3,k), plot(t,Yc(2:(T+1),i,j),'-r')
            hold on
            subplot(ceil(count/3),3,k), plot(t,Yn(2:(T+1),i,j),'--b')
            hold on
            subplot(ceil(count/3),3,k), plot(t,Yc(2:(T+1),i,j) + Yn(2:(T+1),i,j),':k')
            hold off
            %axis([0 T 0 max(max(Ac + An))])
            xlabel('time');
            ylabel(['y(',num2str(i),',',num2str(j),')']);
        end
    end
end
legend('compliant vehicles','non-compliant vehicles','compliant and non-compliant vehicles')
[ax,h3]=suplabel('Flows between cells','t');
set(h3,'FontSize',18)

%plot xn, xc for each cell for all times T
figure
k = 0;
t = 1:T;
for i = 1:numcells
    k = k + 1;
    subplot(ceil(count/3),3,k), plot(t,xc(2:(T+1),i),'-r')
    hold on
    subplot(ceil(count/3),3,k), plot(t,xn(2:(T+1),i),'--b')
    hold on
    subplot(ceil(count/3),3,k), plot(t,xc(2:(T+1),i) + xn(2:(T+1),i),':k')
    hold off
%         if i == 1
%             axis([0 T 0 max(zc + zn)])
%         elseif i == numcells
%             axis([0 T 0 max(zc(numcells:numcells:T*numcells) + zn(numcells:numcells:T*numcells))])
%         else        
%             axis([0 T 0 13])
%         end    
    xlabel('time');
    ylabel(['x(',num2str(i),')']);
end
legend('compliant vehicles','non-compliant vehicles','compliant and non-compliant vehicles')
[ax,h3]=suplabel('Number of vehicles in cell','t');
set(h3,'FontSize',18)

%plot the error (error = LH - RH) for each cell for all times T
figure
k = 0;
t = 1:T;
for i = 1:numcells
    k = k + 1;
    subplot(ceil(count/3),3,k), plot(t,error(k:numcells:T*numcells),'r')
    %axis([0 T 0 max(error)])
    xlabel('time');
    ylabel(['Error for cell ',num2str(i)]);
end
[ax,h3]=suplabel('FIFO violation at cell','t');
set(h3,'FontSize',18)

%==========================================================================
%Testing whether the total flow-in is equal to the total flow-out of the
%network.
[simulation, T_multiplicator] = flowin_flowout_check(t_mult_start, compliant_demand, noncompliant_demand, ...
    compliant_startvehicles, noncompliant_startvehicles, sink, Yn, Yc);

end