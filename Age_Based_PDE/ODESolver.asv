function dYdt = ODESolver(t,Y,n, beta, gamma, u, sigma)
 lOI = 1/180; %% 1/number of days immune- loss of immunity rate
 
S = Y(1:n)';   %% Susceptibles
I = Y(n+1:(2*n))' ;    %infected populations governed by choice of n
R = Y(2*n+1);   %% Recovered
N = sum(S)+ sum(I)+ R; %% Total Population

dSdt = - (S'.*beta*sum(I)); % evolution of susceptible. for seasonal add "+lOI*R" and define loI- loss of immunity rate
dIdt = S'.*beta*sum(I) - gamma .* I' ;% evolution of Infected populations (matrix form)
dRdt = sum(gamma .* I'); % Recovered for seasonal add "-lOI*R"
dYdt = [dSdt ;  dIdt; dRdt];% Solution matrix
