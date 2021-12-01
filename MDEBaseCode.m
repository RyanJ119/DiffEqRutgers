%   Model for progression of disease with variants
clear all;
%parameters: vas (number of variants acting on the system), n (mutations of variants), R0(replications rate), Beta(infection
%rate), Gamma(recovery rate), markov chain (mutation factor),  

vas = 5; %number of variants affecting the system: variants associated with specific data used to build R function
n = 1000; %number of different mutations of variants n>1 (use multiples of vas)
mu = linspace(0,1,n) ; %x range of distributions: taking 0-1 and placing n evenly distributed points 
initI = zeros(1, n); %initialize the infected populations with zeroes

initI(n/2) = 1; %set middle element to one infected
%initI(1) = 1; 
 for i = 1:vas
 initI((i/vas)*n) = 1;     %update the initial infected populations using values or functions
 end

%initI = variants/100+1;% this will need to be fixed to be a function matching the four variants

Yo = [600000; initI' ;0];%% Initial S, I1, I2,... In, R


S = []; %susceptible populations
I = [];% infected pop
R = [];%recovered pop
Iend=[];
% T = [0.9 0.1 0.00 0.00; %A
%     0.02 0.9 0.03 0.05; %B
%     0.02 0.03 0.9 0.05; %G
%     0.02 0.03 0.05 .9]; %D








    Ro= mu/3+1.5 ; % replication rate function 
    gamma = (mu+1)./(mu+1) -13/14; % recovery rate (here it is just 1/14
    beta=Ro.*gamma; %find beta for all infected populations
   
   

    
daysUpdate = 2; % Number of days between mutations (swapping between infected groups) 
totalDays = 500; %total days of program
i = totalDays/daysUpdate;  % number of times to run the solver


 %%%%%%%%%%%%%%%%%% Call the solver i times with a transition of infected inbetween calls
 for c = 1:i
     tRange = daysUpdate*(i-1):1:daysUpdate*i;   %% number of days to run before breaking out of solver to mutate infected
    
%     if  length(R)>300
%         
%         Ro=zeros(n);
%         beta=Ro.*gamma; %%%%Use this to change replication rate mid run
%         
%     end
    

    
    [tSol,YSol] = ode45(@(t,Y) SIRmodels(t,Y, n, beta', gamma'), tRange, Yo);
%%%% concatinating matrices double counts last previous entry/first new
%%%% entry. one of them here, then concatonate old solution with new
%%%% solution
    if length(S)>0 
        S(end)=[];
    end
    
    S = vertcat(S, YSol(:,1));
 
    if length(I)>0
        I(end, :) = [] ;
    end

    I = vertcat(I,  YSol(:,2:(n+1)));

    if length(R)>0
        R (end)=[];
    end
    
    R = vertcat(R, YSol(:,n+2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update step
m = 20;
delv = 1/10;
delm = 1/m;

flip = I';              %Turn I in order to do mutation step
Iend = flip(:, end);
probdist = Iend'./sum(Iend);
F = zeros(1, n);
F2 = [];
v = [];
 for k = 1:n 
        
        if probdist(1,k) > 0
            F(k) = sum(probdist(1:k));                     %should F definitely be cumulative distribution including zeros?
        end                                     % Should F be calculated fully before the update sequence, or during?
    end
    for k = 1:n
        
        if  F(k) > 0
            
            for j = 1:floor(probdist(k)/delm)         %Does not happen when I(k) < .05
                %              if k == 1
                %                 F2(k,j) = 0;
                %              end
                
                if j == floor(probdist(k)/delm) && k>1 && j>1
                    F2(k,j) = F(k-1) + F(k) - F2(k, j-1);
                end
                if  k>1
                    F2(k,j) = F(k-1)+delm*F(k);
                end
                
                if k > 1
                    
                    v(k,j) = floor((1/delv)*rho(F(k-1) + sum(F2(k,1:j))));
                    e1 = zeros(1, n);
                    e2 = zeros(1, n);
                    e1(k) = 1;
                    if probdist(k) -  e1*F2(k,j) >=0
                        
                        if (k+v(k,j))<=n && (k+v(k,j)>=1)
                            e2(k+v(k,j)) = 1;
                            probdist = probdist -  e1*F2(k,j) + e2*F2(k,j);
                            %             else
                            %                 w = 5
                            
                        elseif (k+v(k,j))>n
                            e2(n) = 1;
                            probdist = probdist -  e1*F2(k,j) + e2*F2(k,j);
                        elseif (k+v(k,j))<1
                            e2(1) = 1;
                            probdist = probdist -  e1*F2(k,j) + e2*F2(k,j);
                        end
                    end
                end
            end
        end
    end

    
    
    

Iend = probdist.*sum(Iend);

        



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



%  figure;
%  graphplot(mc,'LabelEdges',true);       %markov chain visual
%  
% figure;
% imagesc(Trmc);
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

sumItotal = 0; %set the total infected to zero before summing
 for i = 1:n
%     plot(tSol,I(:,i), 'yellow');  Plotting all mutations of variants if desired
     sumItotal = sumItotal+I(:,i);  %sum Infected for total infected number
 end

sumI = zeros(n/vas, totalDays+1)';   %set up a matrix to conatenate our infected populations into specific variants


for j = 1:vas
    for i = 1:(n/vas)
        
        sumI(:,j) = sumI(:,j)+I(:,(j-1)*(n/vas)+i);  %concatenate n infected populations into vas categories
       
    end
 
end

figure;
h(1) = plot(tSol,S, 'DisplayName', 'susceptible');
hold on;
for i = 1:vas
h(i+1) =     plot(tSol, sumI(:,i),'DisplayName', ['Infection ' num2str(i)]); %plot new infected;
end

h(vas+2) = plot(tSol, sumItotal, 'DisplayName', 'Total Infected ') ;

h(vas+3) = plot(tSol,R, 'DisplayName', 'Recovered ');

legend(h, 'FontSize', 18)
xlabel("Days", 'FontSize', 18) 
ylabel("Number of Individuals", 'FontSize', 18)
h = [];
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

function [result] = rho(x)

result = 2*x-1;

end