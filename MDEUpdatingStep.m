clear all;
n = 100; %update next piece of mui 
m = 20;
v = .5*m;
delv = 1/v;
delm = 1/m;
mu = zeros(n); %initialize the infected populations with zeroes
%mu(n/4) = .42;
mu(1, n/2) = 1;
%mu(n/2) = .58; %set middle element to one infected

for k = 1:80
    
%     for i = 1:n
%  
%          F(i) = sum(mu(1:i));
%         end 
   
    
    for i = 1:n
        F(i) = sum(mu(k,1:i));
        
        
        
        for j = 1:ceil(mu(k,i)/delm)
            
            if mu(k, i) >=delm
                
                
                %pieces = zeros(1,j);
                % F2 = zeros(1,j);
                if      j ~= ceil(mu(k, i)/delm) && j>1
                    F2(j) = F2(j-1)+delm;
                    pieces(j) = delm;
                elseif      j == 1
                    F2(j) = delm;
                    pieces(j) = delm;
                elseif j == ceil(mu(k, i)/delm)
                    F2(j) = mu(k, i);
                    pieces(j) =mu(k, i)-F2(j-1) ;
                end
            end
            if mu(k, i) <=delm
                F2(j) = mu(k, i);
                pieces(j) = mu(k, i);
            end
     
        if k == 60
        if i == 61
                phi(F(i-1)+F2(j))
        end
        end 
            %vel = zeros(1,j);
            if i > 1
                vel(j) = floor((1/delv)*phi(F(i-1)+F2(j))+0.000001);
                if i == 61
                phi(F(i-1)+F2(j))
                end
            elseif i == 1
                vel(j) = 0;
            end
            
            e1 = zeros(1, n);
            e2 = zeros(1, n);
            e1(i) = 1;
            
            if (i+vel(j))<=n && (i+vel(j)>=1)
                
                e2(i+vel(j)) = 1;
                
                %  e2
            elseif (i+vel(j))>n
                e2(n) = 1;
                
            elseif (i+vel(j))<1
                e2(1) = 1;
                
            end
            mu(k+1,:) =  mu(k+1,:)+e2*pieces(j);
            
            %         if  e2 ~=zeros(1,n)
            %             e2
            %         end
        end
    end
    for i = 1:n
        if mu(k, i) < .00005
            mu(k, i) = 0;
        end
        
    end
end

function [result] = phi(x)

result = 2*x-1;

end

   
