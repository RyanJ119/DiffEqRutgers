%   Model for progression of disease with variants
clear all;
%parameters: vas (number of variants acting on the system), n (mutations of variants), R0(replications rate), Beta(infection
%rate), Gamma(recovery rate), markov chain (mutation factor),

vas = 5; %number of variants affecting the system: variants associated with specific data used to build R function
n = 1000; %number of different mutations of variants n>1 (use multiples of vas)
mu = linspace(0,1,n) ; %x range of distributions: taking 0-1 and placing n evenly distributed points
initI = zeros(1, n); %initialize the infected populations with zeroes
% for i = 1:n

% initI(i) = 1; %set middle element to one infected
% end
initI(n/2) = 1;
%  for i = 1:vas
%  initI((i/vas)*n) = 1;     %update the initial infected populations using values or functions
%  end

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








Ro=3*(mu-.5).^2+1.25 ; % replication rate function
gamma = (mu+1)./(mu+1) -13/14; % recovery rate (here it is just 1/14
beta=Ro.*gamma; %find beta for all infected populations




daysUpdate = 2; % Number of days between mutations (swapping between infected groups)
totalDays = 400; %total days of program
i = totalDays/daysUpdate;  % number of times to run the solver
probdist = zeros(i+1,n);

%%%%%%%%%%%%%%%%%% Call the solver i times 
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
    %%%% entry. Delete one of them here, then concatonate old solution with new
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
    
    
    %%%%%%%%%%%%%%%%%% Call the solver i times 
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update step
   
    m = 20;
    v = .5*m;
    delv = 1/v;
    delm = 1/m;
    
    flip = I';              %Turn I in order to do mutation step
    Iend = flip(:, end);
    
    
    probdist(c,:) = Iend'./sum(Iend);
     if  probdist(c,1) ==0
    F = zeros(n);
    
    for k = 1:n
        F(k) = sum(probdist(c,1:k));
        
        
        
        for j = 1:ceil(probdist(c,k)/delm)
            
            if probdist(c, k) >=delm
                
                
                
                if      j ~= ceil(probdist(c, k)/delm) && j>1
                    F2(j) = F2(j-1)+delm;
                    pieces(j) = delm;
                elseif      j == 1
                    F2(j) = delm;
                    pieces(j) = delm;
                elseif j == ceil(probdist(c, k)/delm)           %Split mass at a point into pieces sized delm 
                    F2(j) = probdist(c, k);
                    pieces(j) =probdist(c, k)-F2(j-1) ;
                end
            end
            if probdist(c, k) <=delm
                F2(j) = probdist(c, k);
                pieces(j) = probdist(c, k);
            end
            
            
            %vel = zeros(1,j);
            if k > 1
                vel(j) = floor((1/delv)*phi(F(k-1)+F2(j))+0.000001);
                
            elseif k == 1
                vel(j) = 0;
            end
            
            e1 = zeros(1, n);
            e2 = zeros(1, n);
            e1(k) = 1;
            
            if (k+vel(j))<=n && (k+vel(j)>=1)
                
                e2(k+vel(j)) = 1;
                
                %  e2
            elseif (k+vel(j))>n
                e2(n) = 1;
                
            elseif (k+vel(j))<1
                e2(1) = 1;
                
            end
            probdist(c+1,:) =  probdist(c+1,:)+e2*pieces(j);
            
            %         if  e2 ~=zeros(1,n)
            %             e2
            %         end
        end
    end
    
    for w = 1:i+1
        for q = 1:n
            if probdist(i+1, n) < .0000005
                probdist(i+1, n) = 0;
            end
        end
    end
    
    
   
    Iend = probdist(c+1,:).*sum(Iend);
    
    
    
    
    end
    I(end, :) = Iend;
    
    Yo = [S(end); I(end, :)' ;R(end)];

end






    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update step






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting

tRange = 0:1:length(S)-1;
tSol = tRange;                  %Define time range  for plotting




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


figure('name','dynamicsOfMDECode'); %Plotting the dynamics altogether

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

%h = [];




figure('name','input values'); %Plotting R Beta and Gamma

subplot(2,2,1);
plot(mu,gamma)
title('Gamma')
subplot(2,2,2);
plot(mu,Ro)
title('Replication Rate')
subplot(2,2,3);
plot(mu,beta)
title('Beta')
subplot(2,2,4); 
plot(mu,phi(mu))
title('Phi')

figure;
 ss = size(I);
 [x,y] = ndgrid(1:ss(2),1:ss(1));

surf(y,x,I.');
colormap turbo
c = colorbar
shading interp
xlabel('Days')
ylabel('Mutations')
zlabel('Infected')

figure;
 ss = size(probdist);
 [x,y] = ndgrid(1:ss(2),1:ss(1));

surf(y,x,probdist.');
colormap turbo
c = colorbar
shading interp
xlabel('Days')
ylabel('Mutations')
zlabel('Percent of total infections')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Functions

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



function [result] = phi(x)

result = 6*(x-.5).*(x-.5).*(x-.5);
%result = -1*sin(10*x)*.4-.1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Functions