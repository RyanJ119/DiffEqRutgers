function dYdt = ODESolver(t,Y,n, beta, gamma, u, sigma)
 lOI = 1/180; %% 1/number of days immune- loss of immunity rate
 
S = Y(1);   %% Susceptibles
I = Y(2:(n+1))' ;    %infected populations governed by choice of n
R = Y(n+2);   %% Recovered
H = Y(n+3);
N = S+ sum(I)+ R+H; %% Total Population

dSdt = -sum(beta * (S/N) .* I')*u; % evolution of susceptible. for seasonal add "+lOI*R" and define loI- loss of immunity rate
dIdt = u*beta * (S/N) .* I' - gamma .* I' - sigma .* I';% evolution of Infected populations (matrix form)
dRdt = sum(gamma .* I') +H/14; % Recovered for seasonal add "-lOI*R"
dHdt = sum(sigma .* I') -H/14 ;
dYdt = [dSdt ;  dIdt; dRdt; dHdt];% Solution matrix
