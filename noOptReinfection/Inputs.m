%   Model for progression of disease with variants

%input parameters: vas (number of variants acting on the system), n (mutations of variants), R0(replications rate), Beta(infection
%rate), Gamma(recovery rate), markov chain (mutation factor),

vas = 10; %number of variants affecting the system: variants associated with specific data used to build R function
n = 1000; %number of different mutations of variants n>1 (use multiples of vas)
mu = linspace(0,1,n) ; %x range of distributions: taking 0-1 and placing n evenly distributed points
initI = zeros(1, n); %initialize the infected populations with zeroes
initSr = zeros(1, n);
initR = zeros(1, n);
% for i = 1:n

% initI(i) = 1; %set middle element to one infected
% end
initI(n/2) = 1;

%  for i = 1:vas
%  initI((i/vas)*n) = 1;     %update the initial infected populations using values or functions
%  end

%initI = variants/100+1;% this will need to be fixed to be a function matching the four variants
initialS = 8.882e6;

Yo = [initialS; initI' ;initR'; initSr' ];%% Initial S, I1, I2,... In, R, H


S = []; %susceptible populations
I = [];% infected pop
R = [];%recovered pop
Sr = [];
Iend=[];
% T = [0.9 0.1 0.00 0.00; %A
%     0.02 0.9 0.03 0.05; %B
%     0.02 0.03 0.9 0.05; %G
%     0.02 0.03 0.05 .9]; %D

Ro=-.3*sin(100*mu/pi-2.1)+1.8 ;
%Ro=(mu+1)./(mu+1) +1;
daysBetweenGovtUpdates = 30;

%Ro=(mu-.5).^2+1.75 ; % replication rate function
gammaRate = (mu+1)./(mu+1) -13/14; % recovery rate (here it is just 1/14

sigma = (mu+1)./(mu+1) -10/11; % 1/30 
sigma = sigma';


betaRate=Ro.*gammaRate; %find betaRate for all infected populations

betaHat = zeros(n, n);
 for i= 1:length(betaRate)
     betaHat(i,:) = betaRate;
     for j= 1: 1:length(betaRate)
         if  i+j<length(betaRate) && abs(i-j)<= 80
             betaHat(i,j) = 0;
         end
    end
 end

betaRate = betaRate';
gammaRate = gammaRate';
%betaHat = betaHat;
%sigma= 1/40; %reinfection rate




daysUpdate = 2; % Number of days between mutations (swapping between infected groups)
totalDays = 800; %total days of program
runSolver = totalDays/daysUpdate;  % number of times to run the solver
probdist = zeros(runSolver+1,n);


%u = ones(1,runSolver);

%u = .99; %1 = no lockdown 0 = full lockdown
%% Mutating constants 

    m = 40;
    v = .1*m;       %Speed at which mutation will travel 
    delv = 1/v;
    delm = 1/m;     %Size of one piece of infected to be mutated
    
