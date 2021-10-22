%   Model for progression of disease with variants

%parameters: n (number of variants), R0(replications rate), Beta(infection
%rate), Gamma(recovery rate), markov chain (mutation factor) 

vas = 2; %number of variants affecting the system: variants associated with specific data used to build R function
n = 16; %number of variants n>1 (use multiples of vas)
variants = linspace(0,1,n) ; %x range of distributions: taking 0-4 
initI = variants/100+1; % this will need to be fixed to be a function matching the four variants

Yo = [600000; initI' ;0];%% Initial S, I1, I2,... In, R
sumItotal = 0;


Trmc = []; % transition matrix
S = []; %susceptible populations
I = [];% infected pop
R = [];%recovered pop
Iend=[];
% T = [0.9 0.1 0.00 0.00; %A
%     0.02 0.9 0.03 0.05; %B
%     0.02 0.03 0.9 0.05; %G
%     0.02 0.03 0.05 .9]; %D




%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = NaN(n);

for i = 1:n
    for j = 1:n
        if i == j
            P(i,i) = 0.9; %choose percent of each infected that stay in their own cateogry each update
            
        elseif abs(i-j) ~= 1
            P(i,j) = 0;     %Only Infected populations next to each other feed to each other
        end
    end
end
                                        % Produce a random transition
                                        % matrix with 0.9 on the diagonal
                                        % and 0 outside of he off diagonal 
                                        % and n states
rng(1); % For reproducibility
mc = mcmix(n,'Fix',P); % add ',Zeros',3 for random zeroes 

Trmc = mc.P;

%%%%%%%%%%%%%%%%%%%%%%%%%



    Ro= variants/3+1.5 ; % replication rate function 
    gamma = (variants+1)./(variants+1) -13/14; % recovery rate (here it is just 1/14
    beta=Ro.*gamma; %find beta for all infected populations
   
   
%  Ro=[1.35, 1.2, 1.6, 2.1];  %Replication  Rates of variants	
 % gamma = [1/14,1/14,1/14,1/14] ;  %% Choose Recovery Rate (Gamma(I_i)) should be 0.072 - 0.142 (1/ 7-14 days)
 % beta=Ro.*gamma; %find beta for all infected populations


    
daysUpdate = 7; % Number of days between mutations (swapping between infected groups) 
totalDays = 600; %total days of program
i = totalDays/daysUpdate;  % number of times to run the solver


 %%%%%%%%%%%%%%%%%% Call the solver i times with a transition of infected inbetween calls
 for c = 1:i
     tRange = daysUpdate*(i-1):1:daysUpdate*i;   %% number of days to run before breaking out of solver to mutate infected
    
%     if  length(R)>300
%         
%         Ro=zeros(n);
%         beta=Ro.*gamma; %find beta for all infected populations
%         
%     end
    

    
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
    R = vertcat(R, YSol(:,n+2));

    flip = I';
    Iend = flip(:, end);


    Iend = ((Iend'./sum(Iend))*Trmc).*sum(Iend); %find probability distribution, multiply by transition matrix,then multiply by total infected again

    I(end, :) = Iend;
    Yo = [S(end); I(end, :)' ;R(end)];
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
tRange = 0:1:length(S)-1;
tSol = tRange;                  %Define time range  for plotting



%%%%%%%%% Plot the Markov Chain Defined Above
namesmark = strings;
names = strings;
for i = 1:n
    namesmark(i) = ['Infection ' num2str(i)]  ; %Name infected populations (lease infective to most)
end 


mc = dtmc(Trmc,'StateNames',[namesmark]); % define markov chain for plotting 

for i = 1:n/vas
    names(i+1) = ['Infection ' num2str(i)]  ;  
end 
                                                % add susceptible, total infected and
                                                % recovered for plotting 
names (1) = "susceptible";
names (n+2) = "total infected";
names (n+3) = "recovered";



%  figure;
%  graphplot(mc,'LabelEdges',true);       %markov chain visual
%  
% figure;
% imagesc(T);
% colormap(jet);
% colorbar;
% axis square
% h = gca;                    %heatmap visual for markov chain
% h.XTick = 1:n;
% h.YTick = 1:n;
% title 'Transition Matrix Heatmap';
% hold on;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting solution to SIR model
figure;
plot(tSol,S)
hold on
for i = 1:n
    plot(tSol,I(:,i), 'yellow');
    sumItotal = sumItotal+I(:,i);
     
end
sumI = zeros(n/vas, length(I))';   %set up a matrix to conatenate our infected populations into specific variants


for j = 1:vas
    for i = 1:(n/vas)
    
    
        sumI(:,j) = sumI(:,j)+I(:,(j-1)*(n/vas)+i);  %concatenate n infected populations into vas categories
       
    end
 
end
 

for i = 1:n/vas
    plot(tSol, sumI(:,i), 'blue') %plot new infected
end

plot(tSol, sumItotal, 'red')                    
plot(tSol,R)
legend(names, 'FontSize', 18)
xlabel("Days", 'FontSize', 18) 
ylabel("Number of Individuals", 'FontSize', 18)

function dYdt = SIRmodels(t,Y,n, beta, gamma)
   % lOI = .02
    S = Y(1);   %% Susceptibles
    I = Y(2:(n+1))' ;    %infected populations governed by choice of n 
    R = Y(n+2);   %% Recovered
    N = S+ sum(I)+ R; %% Total Population
    
    dSdt = -sum(beta * (S/N) .* I'); % evolution of susceptible. for seasonal add "+lOI*R" and define loI- loss of immunity rate
    dIdt = beta * (S/N) .* I' - gamma .* I';% evolution of Infected populations (matrix form)
    dRdt = sum(gamma .* I') ; % Recovered for seasonal add "-lOI*R"
    
    dYdt = [dSdt ;  dIdt; dRdt];% Solution matrix 
end