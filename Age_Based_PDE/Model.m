clear all;

Inputs; %bring in all chosen input values

%% Call the solver i times 
for c = 1:i
    
    tRange = daysUpdate*(i-1):1:daysUpdate*i;   %% number of days to run before breaking out of solver to mutate infected
    
    %     if  length(R)>300
    %
    %         Ro=zeros(n);
    %         beta=Ro.*gamma; %%%%Use this to change replication rate mid run
    %
    %     end
    
    
    
    [tSol,YSol] = ode45(@(t,Y) ODESolver(t,Y, n, beta', gamma', u, sigma'), tRange, Yo);
  
    %%%% concatinating matrices double counts last previous entry/first new
    %%%% entry. Delete one of them here, then concatonate old solution with new
    %%%% solution
    
    if length(S)>0
        S(end,:)=[];
    end
    
    S = vertcat(S, YSol(:,1:n));
    
    if length(I)>0
        I(end, :) = [] ;
    end
    
    I = vertcat(I,  YSol(:,n+1:2*n));
    
    if length(R)>0
        R (end)=[];
    end
    R = vertcat(R, YSol(:,(2*n+1)));
    
    

   
end






Plotting; %plot all returned values
