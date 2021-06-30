% Model parameters

% to test for a good choice of beta and gamma, beta/gamma should be about
% 1-5


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Define Markov Chain for Infected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Interactions
Tr = [0.9 0.05 0.02 0.03; % Mutations probability for infected1
     0.02 0.9 0.03 0.05;% Mutations probability for infected2
     0.02 0.03 0.9 0.05;% Mutations probability for infected3
     0.02 0.03 0.05 0.9]; % Mutations probability for infected4
 
 
 %Transition Matrix Print
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;

imagesc(Tr);
colormap(jet);
colorbar;
axis square
h = gca;
h.XTick = 1:4;
h.YTick = 1:4;
title 'Transition Matrix Heatmap';
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta = 0.0; % rate of immunity loss
N = 60000000; % Total population N = S + I + R


I10 = 100; % initial number of infected
I20 = 1; % initial number of infected
I30 = 1; % initial number of infected
I40 = 1; % initial number of infected


T = 180; % period of 180 days
dt = 1; % time interval for updating

% Calculate the model
[S,I1, I2, I3, I4,R, Tr] = sir_model(delta,N,I10, I20, I30, I40,Tr, T,dt);
% Plots that display the epidemic outbreak
%fprintf('Value of parameter R0 is %.2f',N*beta/gamma)
t = 0:dt:T-dt;
% Curve

figure;
plot(t,S,'b',t,I1,t,I2,t,I3,t,I4,'r',t,R,'g','LineWidth',2); hold on;
xlabel('Days'); ylabel('Number of individuals');
legend('S','I1','I2','I3','I4','R');
hold on;
% Map
function [S,I1, I2, I3, I4,R,Tr, T] = sir_model(delta,N,I10,I20,I30,I40,Tr, T,dt)
    % if delta = 0 we assume a model without immunity loss
    
    
    
    S = zeros(1,T/dt);
    S(1) = N;
    I1 = zeros(1,T/dt);
    I1(1) = I10;
    I2 = zeros(1,T/dt);
    I2(1) = I20;
    I3 = zeros(1,T/dt);
    I3(1) = I30;
    I4 = zeros(1,T/dt);
    I4(1) = I40;
    R = zeros(1,T/dt);

    for t = 1:(T/dt)-1
        
% to test for a good choice of beta and gamma, beta/gamma should be about
% 1-5
%0.3+0.5*sin(t)
%

        beta1 = .3; 
        beta2 = .3; 
        beta3 = .3; %% Infection Rates (Beta(I_i))
        beta4 = .3; 
    
        gamma1 = 0.12;
        gamma2 = 0.12; 
        gamma3 = 0.12;  %% Recovery Rate (Gamma(I_i))
        gamma4 = 0.12; 
        
        gamma = [gamma1 gamma2 gamma3 gamma4];  %vectorize recovery rates
        beta = [beta1 beta2 beta3 beta4];%vectorize infection rates
        I = [I1(t) I2(t) I3(t) I4(t)];

        if t~=0
            if mod(t,50) == 0
                
                I = ((I/sum(I))*Tr)*sum(I); %find probability distribution, multiply by transition matrix,then multiply by total infected again 
                I1(t) = I(1);
                I2(t) = I(2);
                I3(t) = I(3);
                I4(t) = I(4);
            end
        end
        % Equations of the model
        dS  = (-sum(beta * (S(t)/N).* I))*dt; % evolution of susceptible 
        %dS = (-beta1*I1(t)*S(t)) * dt
        dI1 = (beta1*I1(t)*S(t)/N - gamma1*I1(t)) * dt;
        dI2 = (beta2*I2(t)*S(t)/N - gamma2*I2(t)) * dt;
        dI3 = (beta3*I3(t)*S(t)/N - gamma3*I3(t)) * dt;
        dI4 = (beta4*I4(t)*S(t)/N - gamma4*I4(t)) * dt;
        dR =  sum(gamma .* I)* dt; % Recovered  * dt;
        
        S(t+1) = S(t) + dS;
        I1(t+1) = I1(t) + dI1;
        I2(t+1) = I2(t) + dI2;
        I3(t+1) = I3(t) + dI3;
        I4(t+1) = I4(t) + dI4;
        R(t+1) = R(t) + dR;
    end
end

