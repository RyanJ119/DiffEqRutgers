
function dYdt = ODESolver(t,Y,n, beta, gamma, betaHat, sigma)
 %lOI = 1/180; %% 1/number of days immune- loss of immunity rate
 
S = Y(1);   %% Susceptibles
I = Y(2:(n+1))  ;    %infected populations governed by choice of n

R = Y( (n+2) : 2*n+1 );   %% Recovered

Sr = Y( 2*n+2 : 3*n+1 ); 

N = S+ sum(I)+ sum(R)+sum(Sr); %% Total Population

dSdt = -sum(beta * (S/N) .* I); % evolution of susceptible. for seasonal add "+lOI*R" and define loI- loss of immunity rate
dIdt = beta * (S/N) .*  I + sum((betaHat.*Sr.*I), 2)./N - gamma .*  I   ;% evolution of Infected populations (matrix form)
dRdt = gamma .*  I - sigma.*R; % Recovered for seasonal add "-lOI*R"
Srdt =  sigma.*R - sum((betaHat.*Sr.* I), 2)./N ;   % mtimes( betaHat, Sr' )/N .* I'

dYdt = [dSdt ;  dIdt; dRdt; Srdt ];% Solution matrix

