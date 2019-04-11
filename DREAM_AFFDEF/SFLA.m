clear all
clc
disp ('created by sadegh sadeghitabas')
tic
global ss PARSout NSout Spar Qobs num sst PARt NSt
num=24;
ss=1;
sst=1;
%shuffled frog-leaping algorithm
%initialize the five parameters of SFLA
m=3; %the number of memeplexes
n=3; %the number of frogs in a memeplex
q=2; %the number of frogs in a submemeplex
Ne=5; %the number of evolution steps in a memeplex between two successive shufflings
Smax = 1;  %the maximum step size

%initialize the other parameters
F=m*n;   %the virtual population size,and F=m*n
d=9; %the number of decision variables memotypes in a meme carried by a frog
Pmin =[0.3 300 0.05 0.05 0.5 0.001 20000 0.05 0.05]; %the frogs' d-dimensions maximum position
Pmax =[1.5 100000 10 10 20 0.1 800000 0.9 0.7];%the frogs' d-dimensions minimum position
%generation number
MAXGEN=5;

% generate the vitual population
for i=1:F
    U(i,:)=Pmin+(rand(1,d).*(Pmax-Pmin));
end

loop = 1;
while loop<=MAXGEN
    
    %compute the performance value of each frog in population
    for i=1:F
        fitness(i)=NS(U(i,:));
    end
    
    %sort the F frogs in order of decreasing performance value
    [fitnessSort,index]=sort(fitness,'ascend');
    for i=1:F
        X(i,:)=U(index(i),:);
    end
    %the best frog's position Px in the entire population
    Px=X(1,:);
    
    %partition frogs into m memeplexes Y1,Y2,...,Ym
    for i=1:m
        for j=1:n
            Y(i,j,:)=X(i+m*(j-1),:);
        end
    end
    %memetic evolution within each memeplex
    for im=1:m
        for iN=1:Ne
            %record the best frog's position as PB
            pb=Y(im,1,:);
            pw=Y(im,n,:);
            %compute the step size
            for i=1:d
                S1(i) = rand*(pb(i)-pw(i));
                if S1(i)>=0
                    S1(i)=min([S1(i),Smax]);
                else
                    S1(i)=max([S1(i),-Smax]);
                end
            end
            for i=1:d
                Zq(i)= pw(i)+S1(i);
            end
            %improve the worst frog's position
            
            
            %if Zq is within the feasible space,compute the new performance
            %value,else chose the step size again
            
            fail=0;
            for i=1:d
                if Zq(i)>Pmax(i)
                    fail=1
                end
                if Zq(i)<Pmin(i)
                    fail=1
                end
            end
            
            
            %Zq is within the feasible space,compute the new performance value
            if fail==0
                Zqfitness = NS(Zq);
                if Zqfitness<NS(Y(im,q,:));
                    Z(im,n,:)=Zq;
                else
                    fail=1;
                end
            end
            
            %if Zq isn't within feasible space or the new performance value
            %isn't better than the old.
            if fail==1
                %compute the step size again
                for i=1:d
                    S2(i) = rand*(Px(i)-pw(i));
                    if S2(i)>=0
                        S2(i)=min([S2(i),Smax]);
                    else
                        S2(i)=max([S2(i),-Smax]);
                    end
                end
                
                %compute the new position
                for i=1:d
                    Zq(i)= pw(i)+S2(i);
                end
                %judge if Zq is within the feasible space
                fail=0;
                for i=1:d
                    
                    if Zq(i)>Pmax(i)
                        fail=1;
                    else
                        if Zq(i)<Pmin(i)
                            fail=1;
                        end
                    end
                end
                
                %if Zq is within the feasible space,compute the new performance
                %value
                if fail==0
                    Zqfitness = NS(Zq);
                    if Zqfitness<NS(Y(im,q,:))
                        Z(im,n,:)=Zq;
                    else
                        fail=1;
                    end
                end
            end
            %if Zq is either infeasible or not better than the old, then
            %generate a new frog randomly at a feasible location
            if fail==1
                Zq=Pmin+(rand(1,d).*(Pmax-Pmin));
                Z(im,n,:)=Zq;
            end
            
            
            %upgrade the memeplex,replace Z in their original location in Yim
            
            Y(im,n,:)=Z(im,n,:);
            
            
            %sort Yim in order of decreasing performance value.
            for i=1:n
                Yfitness(i)=NS(Y(im,i,:));
            end
            [Yfitness,index]=sort(Yfitness,'ascend');
            for i=1:n
                Yim(i,:)=Y(im,index(i),:);
            end
            
            Y(im,:,:)=Yim;
        end
    end
    %shuffle memeplexes
    %replace Y1,Y2,..Ym into U,i.e fresh the total population
    for i=1:m
        for j=1:n
            U(i*j,:)=Y(i,j,:);
        end
    end
    loop=loop+1
        
    bestF(loop)=NS(Px)
    bestPar(loop,:)=Px
end
Function=NS(Px)
parameters=Px
disp ('created by sadegh sadeghitabas')
Ntot = sst-1 
ModPred=AFFDEF(Px);
% Compute the RMSE
RMSE = sqrt ( sum ( ( ModPred - Qobs).^2) / num);

% Now create the parameter and total uncertainty
for j = 1:Ntot
    Smod(:,j) = Spar(:,j) + normrnd(0,RMSE,num,1);
end;



% Now sort to get 95% ranges
for j = 1:num
    a = sort(Spar(j,1:Ntot));
    % And take the 2.5 and 97.5 percentiles
    minpar95(j,1) = a(floor(0.025 * Ntot)); maxpar95(j,1) = a(floor(0.975 * Ntot));
    % Same with total uncertainty
    a = sort(Smod(j,1:Ntot));
    % And take the 2.5 and 97.5 percentiles
    mintot95(j,1) = a(floor(0.025 * Ntot)); maxtot95(j,1) = a(floor(0.975 * Ntot));
end;

% Now plot the results 
plot([1:num],minpar95,[1:num],maxpar95); % These are the upper and lower parameter uncertainty bounds
hold on; plot([1:num],mintot95,'R',[1:num],maxtot95,'G'); % These are the upper and lower total uncertainty bounds
plot([1:num],Qobs,'r+');

% Now calculate percentage inside the 95% bound --> should be about 95%
Contained = 1 - length(find(Qobs<mintot95 |  Qobs> maxtot95))/num
Contained1 = 1 - length(find(Qobs<minpar95 |  Qobs> maxpar95))/num
toc
% Note the ranges do not look very good --> because of homoscedasticity. If
% the error is assumed to increase with flow level results look much better