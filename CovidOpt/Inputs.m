%   Model for progression of disease with variants

%input parameters: vas (number of variants acting on the system), n (mutations of variants), R0(replications rate), Beta(infection
%rate), Gamma(recovery rate), markov chain (mutation factor),

vas = 5; %number of variants affecting the system: variants associated with specific data used to build R function
n = 2000; %number of different mutations of variants n>1 (use multiples of vas)
mu = linspace(0,1,n) ; %x range of distributions: taking 0-1 and placing n evenly distributed points
initI = zeros(1, n); %initialize the infected populations with zeroes
% for i = 1:n
u = .8; %1 = no lockdown 0 = full lockdown
% initI(i) = 1; %set middle element to one infected
% end
initI(n/2) = 1;
%  for i = 1:vas
%  initI((i/vas)*n) = 1;     %update the initial infected populations using values or functions
%  end

%initI = variants/100+1;% this will need to be fixed to be a function matching the four variants

Yo = [329500000; initI' ;0; 0];%% Initial S, I1, I2,... In, R, H


S = []; %susceptible populations
I = [];% infected pop
R = [];%recovered pop
H = []; %hospitalized pop
Iend=[];
% T = [0.9 0.1 0.00 0.00; %A
%     0.02 0.9 0.03 0.05; %B
%     0.02 0.03 0.9 0.05; %G
%     0.02 0.03 0.05 .9]; %D

Ro=-.2*sin(9*mu+23)+2.5 ;

%Ro=(mu-.5).^2+1.75 ; % replication rate function
gamma = (mu+1)./(mu+1) -13/14; % recovery rate (here it is just 1/14
beta=Ro.*gamma; %find beta for all infected populations
sigma= (mu+1)./(mu+1)-1 + .01; %hospitalization rate
%sigma= 0; %hospitalization rate




daysUpdate = 2; % Number of days between mutations (swapping between infected groups)
totalDays = 800; %total days of program
i = totalDays/daysUpdate;  % number of times to run the solver
probdist = zeros(i+1,n);

%% Mutating constants 

    m = 20;
    v = .5*m;       %Speed at which mutation will travel 
    delv = 1/v;
    delm = 1/m;     %Size of one piece of infected to be mutated
    
