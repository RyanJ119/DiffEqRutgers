tRange = [0 180];                               %% Time Range
Yo = [999; 1; 2; 3; 10 ;0];

%% Initial S, I1, I2,... In, R

 
S = YSol(:,1);
muI1 = YSol(:,2);
muI2 = YSol(:,3);
muI3 = YSol(:,4); % Extracting Soutions
muI4 = YSol(:,5);

R = YSol(:,6);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Define Markov Chain for Infected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Interactions
T = [0.9 0.05 0.02 0.03; % Mutations probability for infected1
     0.02 0.9 0.03 0.05;% Mutations probability for infected2
     0.02 0.03 0.9 0.05;% Mutations probability for infected3
     0.02 0.03 0.05 0.9]; % Mutations probability for infected
 %Transition Matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Visualizing transition matrix 
 
% stateNames = ["Regime 1" "Regime 2" "Regime 3" "Regime 4"];
% mc = dtmc(T,'StateNames',stateNames);
[tSol,YSol] = ode45(@SIRmodels, tRange, Yo);
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


%%%%%%%%%%%%%%%%%%%%%%%
% Plotting solution to SIR model
figure;
plot(tSol,S)
hold on
plot(tSol,muI1)
plot(tSol,muI2)
plot(tSol,muI3)
plot(tSol,muI4)
plot(tSol,R)
%  hold off
legend("Susceptible","Infected1", "Infected2","Infected3", "Infected4" ,"Recovered")
xlabel("Days")
ylabel("Number of Individuals")


function dYdt = SIRmodels(t,Y)

    S = Y(1);   %% Susceptibles
    I1 = Y(2);%% Infected1
    I2 = Y(3);%% Infected2
    I3 = Y(4);%% Infected3
    I4 = Y(5);%% Infected4
    R = Y(6);   %% Recovered
T = [0.9 0.05 0.02 0.03; % Mutations probability for infected1
     0.02 0.9 0.03 0.05;% Mutations probability for infected2
     0.02 0.03 0.9 0.05;% Mutations probability for infected3
     0.02 0.03 0.05 0.90]; % Mutations probability for infected

   
    N = S+ I1+ I2+ I3+ I4+ R; %% Total Population
    
    beta1 = 1+0.5*sin(t); 
    beta2 = 1+0.5*sin(t); 
    beta3 = 1+0.5*sin(t); %% Infection Rates (Beta(I_i))
    beta4 = 1+0.5*sin(t); 
    
    gamma1 = 0.005*t;
    gamma2 = 0.005*t; 
    gamma3 = 0.005*t;  %% Recovery Rate (Gamma(I_i))
    gamma4 = 0.005*t; 
   
    muI = [I1 I2 I3 I4]; %  infected in each mutation 
    
if R ~= 0 && mod(t,5) < 1
muI = ((muI/sum(muI))*T)*sum(muI); %find probability distribution, multiply by transition matrix,then multiply by total infected again 
end
    
    gamma = [gamma1 gamma2 gamma3 gamma4];  %vectorize recovery rates
    
    beta = [beta1 beta2 beta3 beta4];%vectorize infection rates
    
    dSdt = -sum(beta * (S/N) .* muI); % evolution of susceptible 
    
    dI1dt = beta1 * (S/N) * muI(1) - gamma1 * muI(1);% evolution of Infected population 1 
    dI2dt = beta2 * (S/N) * muI(2) - gamma2 * muI(2);% evolution of Infected population 2 
    dI3dt = beta3 * (S/N) * muI(3) - gamma3 * muI(3);% evolution of Infected population 3 
    dI4dt = beta4 * (S/N) * muI(4) - gamma4 * muI(4);% evolution of Infected population 4
    dRdt = sum(gamma .* muI); % Recovered 
    
     %muI = ((muI/sum(muI))*T)*sum(muI); %find probability distribution, multiply by transition matrix,then multiply by total infected again 

    dYdt = [dSdt ;  dI1dt; dI2dt; dI3dt; dI4dt; dRdt;];% Solution matrix 
end