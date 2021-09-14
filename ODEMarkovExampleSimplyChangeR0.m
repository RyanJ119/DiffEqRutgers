
Yo = [60000000; 1; 0; 0; 0 ;0];%% Initial S, I1, I2,... In, R

S = [];
muI1 = [];
muI2 = []; 
muI3 = []; % Set up solution matrices 
muI4 = [];
R = [];
 muI = [muI1 muI2 muI3 muI4];

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Define Markov Chain for Infected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Interactions
% T = [0 1 0.0 0.0; % Mutations probability for infected1
%      0.0 0.0 1 0.0;% Mutations probability for infected2
%      0.0 0.0 0.0 1;% Mutations probability for infected3
%      1 0.0 0.0 0.0]; % Mutations probability for infected4
%  %Transition Matrix
T = [0.9 0.02 0.05 0.05;
    0.02 0.9 0.03 0.05;
    0.02 0.03 0.9 0.05;
    0.02 0.03 0.05 .9];

 daysUpdate = 7; % Number of days between mutations (updates) 
 totalDays = 400; %total days of program
 i = totalDays/daysUpdate;  % number of times to run the solver
 %%%%%%%%%%%%%%%%%% Call the solver i times with a transition of infected
 %%%%%%%%%%%%%%%%%% inbetween calls
 for c = 1:i
     
     tRange = daysUpdate*(i-1):1:daysUpdate*i;   %% number of days to run before breaking out of solver to mutate infected
     
   % [tSol,YSol] = ode45(@SIRmodels, tRange, Yo); % call solver if Beta and
   % Gamma are called inside function
    
   %%%%%%%%%%%% Test values 
           % beta1 = .1+c*daysUpdate/(4*totalDays) ; %R0 linear from 1.5-5 spread over number of days
        %beta1 = .35-c*daysUpdate/(4*totalDays); %% linear from 3.8-1.5
       % beta2 = .1+(-1*((daysUpdate*c+totalDays)*(daysUpdate*c-totalDays)))/(totalDays*totalDays*6); %parabolic starting at 2.7 going to 1 at the last day 
        %beta3 = .07; %% linear from 5-1.5
        %beta4 = .1+(-1*((daysUpdate*c+totalDays)*(daysUpdate*c-totalDays)))/(totalDays*totalDays*8);%parabolic starting at 3.5 going to 1.5 at the last day 
     
  %      lOM = 0.03; %loss of immunity rate (for seasonal disease) 

           R1=1.2
           R2=1.1
           R3=1   % Replication Rates
           R4=2.5

        gamma1 = 1/14;
        gamma2 = 1/14; 
        gamma3 = 1/14;  %% Choose Recovery Rate (Gamma(I_i)) should be 0.072 - 0.142 (1/ 7-14 days)
        gamma4 = 1/14; 
        
        beta1 = R1*gamma1;
        beta2 = R2*gamma2; 
        beta3 = R3*gamma3; 
        beta4 = R4*gamma4;
    
    gamma = [gamma1 gamma2 gamma3 gamma4];  %vectorize recovery rates
    beta = [beta1 beta2 beta3 beta4];%vectorize infection rates
    
    
    [tSol,YSol] = ode45(@(t,Y) SIRmodels(t,Y,beta, gamma), tRange, Yo);
    
    if length(S)>0 
     S(end)=[];
    end
 S = vertcat(S, YSol(:,1));
 if length(muI1)>0
  muI1 (end)=[];
  muI2 (end)=[];%%%%%%%%%%%%%%%%%%%%% concatinating matrices double count last previous entry/first new entry. delete one of them here
  muI3 (end)=[];
  muI4 (end)=[];
 end
muI1 = vertcat(muI1, YSol(:,2));
muI2 = vertcat(muI2, YSol(:,3)); 
muI3 = vertcat(muI3, YSol(:,4)); % Extracting Soutions by concatinating new days on to old days 
muI4 = vertcat(muI4, YSol(:,5));
if length(R)>0
  R (end)=[];%%%% concatinating matrices double count last previous entry/first new entry. delete one of them here

 end
R = vertcat(R, YSol(:,6));

  I = [muI1(end) muI2(end) muI3(end) muI4(end)];
I = ((I/sum(I))*T)*sum(I); %find probability distribution, multiply by transition matrix,then multiply by total infected again
muI1(end) = I(1);
muI2(end) = I(2); %%%%%%% redistribute new infected to the end of our solution
muI3(end) = I(3);
muI4(end)= I(4);

Yo = [S(end); muI1(end); muI2(end); muI3(end); muI4(end) ;R(end)];%%%% Set new initial conditions for next run 
    
    
 end

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Visualizing transition matrix 
 
% stateNames = ["Regime 1" "Regime 2" "Regime 3" "Regime 4"];
% mc = dtmc(T,'StateNames',stateNames);
tRange = 0:1:length(S)-1;
tSol = tRange;

%%%%%%%%% Plot the Markov Chain Defined Above

mc = dtmc(T,'StateNames',["Infection 1" "Infection 2" "Infection 3" "Infection 4"])

 figure;
 graphplot(mc,'LabelEdges',true);
figure;


imagesc(T);
colormap(jet);
colorbar;
axis square
h = gca;
h.XTick = 1:4;
h.YTick = 1:4;
title 'Transition Matrix Heatmap';
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
% Plotting solution to SIR model
figure;
plot(tSol,S)
hold on
plot(tSol,muI1)
plot(tSol,muI2)
plot(tSol,muI3)
plot(tSol,muI4)
plot(tSol, muI1+muI2+muI3+muI4)
plot(tSol,R)
%  hold off
legend("Susceptible","Infected1", "Infected2","Infected3", "Infected4", "Total Infected","Recovered", 'FontSize', 18)
xlabel("Days", 'FontSize', 18)
ylabel("Number of Individuals", 'FontSize', 18)


function dYdt = SIRmodels(t,Y, beta, gamma)

    S = Y(1);   %% Susceptibles
    muI1 = Y(2);%% Infected1
    muI2 = Y(3);%% Infected2
    muI3 = Y(4);%% Infected3
    muI4 = Y(5);%% Infected4
    R = Y(6);   %% Recovered
    


   
    N = S+ muI1+ muI2+ muI3+ muI4+ R; %% Total Population
    
 
        beta1 = beta(1); 
        beta2 = beta(2); 
        beta3 = beta(3); %% Choose Infection Rates (Beta(I_i)) Here
        beta4 = beta(4);
%        
       % lOI = 0.003; %loss of immunity rate (for seasonal disease) 

        gamma1 = gamma(1);
        gamma2 = gamma(2); 
        gamma3 = gamma(3);  %% Choose Recovery Rate (Gamma(I_i))
        gamma4 = gamma(4); 
    
   % gamma = [gamma1 gamma2 gamma3 gamma4];  %vectorize recovery rates
    %beta = [beta1 beta2 beta3 beta4];%vectorize infection rates
   
   I = [muI1 muI2 muI3 muI4]; %  infected in each mutation 
 

    
    dSdt = -sum(beta * (S/N) .* I); % evolution of susceptible. for seasonal add "+lOI*R"
    
    dmuI1dt = beta(1) * (S/N) * I(1) - gamma1 * I(1);% evolution of Infected population 1 
    dmuI2dt = beta2 * (S/N) * I(2) - gamma2 * I(2);% evolution of Infected population 2 
    dmuI3dt = beta3 * (S/N) * I(3) - gamma3 * I(3);% evolution of Infected population 3 
    dmuI4dt = beta4 * (S/N) * I(4) - gamma4 * I(4);% evolution of Infected population 4
    
    dRdt = sum(gamma .* I) ; % Recovered for seasonal add "-lOI*R"
    
   % dmuIdt = [dmuI1dt dmuI2dt dmuI3dt dmuI4dt];
    
     %I = ((I/sum(I))*T)*sum(I); %find probability distribution, multiply by transition matrix,then multiply by total infected again 

    dYdt = [dSdt ;  dmuI1dt; dmuI2dt; dmuI3dt; dmuI4dt; dRdt];% Solution matrix 
end
