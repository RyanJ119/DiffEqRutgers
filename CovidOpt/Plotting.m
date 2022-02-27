

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
h(vas+4) = plot(tSol,H, 'DisplayName', 'Hospitalized ');

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

subplot(2,1,1);
 ss = size(I);
 [x,y] = ndgrid(1:ss(2),1:ss(1));

surf(y,x,I.');
colormap turbo
c = colorbar;
%caxis([0 1e4])
shading interp
xlabel('Days')
ylabel('Mutations')
zlabel('Infected')



subplot(2,1,2);
 ss = size(probdist);
 [x,y] = ndgrid(1:ss(2),1:ss(1));

surf(y,x,probdist.');
colormap turbo
c = colorbar;
caxis([0 .1])
shading interp
xlabel('Days')
ylabel('Mutations')
zlabel('Percent of total infections')









