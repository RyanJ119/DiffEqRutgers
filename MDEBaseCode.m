n = 4; %number of variants
variant = linspace(0,4,n) ; %x range of distributions: taking 0-4 
inf = variant/100+1; % this will need to be fixed to be a function matching the four variants

Yo = [600000; inf' ;0];%% Initial S, I1, I2,... In, R
add = 0;
S = [];
I = [];
R = [];
Iend=[];
T = [0.9 0.0 0.08 0.02; %A
    0.02 0.9 0.03 0.05; %B
    0.02 0.03 0.9 0.05; %G
    0.02 0.03 0.05 .9]; %D
daysUpdate = 7; % Number of days between mutations (updates) 
totalDays = 400; %total days of program
i = totalDays/daysUpdate;  % number of times to run the solver
 %%%%%%%%%%%%%%%%%% Call the solver i times with a transition of infected inbetween calls
 for c = 1:i
     tRange = daysUpdate*(i-1):1:daysUpdate*i;   %% number of days to run before breaking out of solver to mutate infected
     Ro=[1.35, 1.2, 1.6, 2.1];  %Replication  Rates of variants	
    if  length(R)>300
        Ro=[0, 0, 0, 0];
        beta=Ro.*gamma; %find beta for all infected populations
    end
    gamma = [1/14,1/14,1/14,1/14] ;  %% Choose Recovery Rate (Gamma(I_i)) should be 0.072 - 0.142 (1/ 7-14 days)
    beta=Ro.*gamma; %find beta for all infected populations
    
    [tSol,YSol] = ode45(@(t,Y) SIRmodels(t,Y, n, beta', gamma'), tRange, Yo);

    if length(S)>0 
        S(end)=[];
    end
    
    S = vertcat(S, YSol(:,1));
 
    if length(I)>0
        I(end, :) = [] ;
    end

    I = vertcat(I,  YSol(:,2:(n+1)));

    if length(R)>0
        R (end)=[];%%%% concatinating matrices double count last previous entry/first new entry. delete one of them here
    end
    R = vertcat(R, YSol(:,6));

    flip = I';
    Iend = flip(:, end);


    Iend = ((Iend'./sum(Iend))*T).*sum(Iend); %find probability distribution, multiply by transition matrix,then multiply by total infected again


    I(end, :) = Iend;
    Yo = [S(end); I(end, :)' ;R(end)];
 end
tRange = 0:1:length(S)-1;
tSol = tRange;
%%%%%%%%% Plot the Markov Chain Defined Above
mc = dtmc(T,'StateNames',["Infection 1" "Infection 2" "Infection 3" "Infection 4"]);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting solution to SIR model
figure;
plot(tSol,S)
hold on
for i = 1:n
    plot(tSol,I(:,i));
    add = add+I(:,i);
end
plot(tSol, add)
plot(tSol,R)
legend("Susceptible","Infected1(Alpha)", "Infected2(Beta)","Infected3(Gamma)", "Infected4(Delta)", "Total Infected","Recovered", 'FontSize', 18)
xlabel("Days", 'FontSize', 18) 
ylabel("Number of Individuals", 'FontSize', 18)

function dYdt = SIRmodels(t,Y,n, beta, gamma)
    S = Y(1);   %% Susceptibles
    I = Y(2:(n+1))' ;  
    R = Y(6);   %% Recovered
    N = S+ sum(I)+ R; %% Total Population
    dSdt = -sum(beta * (S/N) .* I'); % evolution of susceptible. for seasonal add "+lOI*R"
    dIdt = beta * (S/N) .* I' - gamma .* I';% evolution of Infected population 1 
    dRdt = sum(gamma .* I') ; % Recovered for seasonal add "-lOI*R"
    dYdt = [dSdt ;  dIdt; dRdt];% Solution matrix 
end