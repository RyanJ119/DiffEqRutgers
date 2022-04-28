
vas = 10; %number of variants affecting the system: variants associated with specific data used to build R function
tRange = 0:1:length(S)-1;
tSol = tRange;                  %Define time range  for plotting


sumItotal = 0; %set the total infected to zero before summing
for i = 1:n
    %     plot(tSol,I(:,i), 'yellow');  Plotting all mutations of variants if desired
    sumItotal = sumItotal+I(:,i);  %sum Infected for total infected number
end

sumRtotal = 0; %set the total infected to zero before summing
for i = 1:n
    %     plot(tSol,I(:,i), 'yellow');  Plotting all mutations of variants if desired
    sumRtotal = sumRtotal+R(:,i);  %sum Infected for total infected number
end

sumSrtotal = 0; %set the total infected to zero before summing
for i = 1:n
    %     plot(tSol,I(:,i), 'yellow');  Plotting all mutations of variants if desired
    sumSrtotal = sumSrtotal+Sr(:,i);  %sum Infected for total infected number
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

h(vas+3) = plot(tSol,sumRtotal, 'DisplayName', 'Recovered ');
h(vas+4) = plot(tSol,sumSrtotal, 'DisplayName', 'Susceptible Recovered ');

legend(h, 'FontSize', 18)
xlabel("Days", 'FontSize', 18)
ylabel("Number of Individuals", 'FontSize', 18)

%h = [];

figure('name','Input Variables');

subplot(2,2,1);
plot(mu,gammaRate)
title('\gamma(\alpha)')
subplot(2,2,2);
plot(mu,Ro)
title('Replication Rate')
subplot(2,2,3);
plot(mu,betaRate)
title('\beta(\alpha)')
subplot(2,2,4); 
plot(mu,phi(mu))
title('\Phi(\alpha)')

figure('name','Infected Over Time');


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


figure('name','Susceptible Recovered Over Time');


 ss = size(I);
 [x,y] = ndgrid(1:ss(2),1:ss(1));

surf(y,x,Sr.');
colormap turbo
c = colorbar;
%caxis([0 1e4])
shading interp
xlabel('Days')
ylabel('Mutations')
zlabel('Susceptible but Recovered')


figure('name','Recovered Over Time');


 ss = size(I);
 [x,y] = ndgrid(1:ss(2),1:ss(1));

surf(y,x,R.');
colormap turbo
c = colorbar;
caxis([0 1e4])
shading interp
xlabel('Days')
ylabel('Mutations')
zlabel('Recovered')




figure('name','Normalized Infected Over Time');
 ss = size(probdist);
 [x,y] = ndgrid(1:ss(2),1:ss(1));
 
 probmat = repelem(probdist,2, 1);
 
 imagesc(probmat');
myColorMap = jet(256);
%myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar()
caxis([0 .0005])
axis square
h = gca;
% h.XTick = 1:4;
% h.YTick = 1:4;
title 'Probability Distribution';
hold on;

%surf(2*y,x,probdist.');
% colormap turbo
% c = colorbar;
%caxis([0 .1])
% shading interp
% xlabel('Days')
% ylabel('Mutations')
% zlabel('Percent of total infections')









