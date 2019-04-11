clc;
clear all;
tic
npop=20;
nvar=9;
global ss PARSout NSout Spar Qobs num sst NSt PARt
num=24;
ss=1;
sst=1;
maxit=25;
w=1;
wdamp=0.99;

c1=2;
c2=2;

xmin=[0.3 300 0.05 0.05 0.5 0.001 20000 0.05 0.05];
xmax=[1.5 100000 10 10 20 0.1 800000 0.9 0.7];
dx=xmax-xmin;

vmax=0.1*dx;

empty_particle.position=[];
empty_particle.velocity=[];
empty_particle.NS=[];
empty_particle.pbest=[];
empty_particle.pbestNS=[];

particle=repmat(empty_particle,npop,1);

gbest=zeros(maxit,nvar);
gbestNS=zeros(maxit,1);

for it=1:maxit
    if it==1
        gbestNS(1)=inf;
        for i=1:npop
            particle(i).velocity=zeros(1,nvar);
            particle(i).position=xmin+(xmax-xmin).*rand(1,nvar);
            particle(i).NS=NS(particle(i).position);
            
            particle(i).pbest=particle(i).position;
            particle(i).pbestNS=particle(i).NS;
            
            if particle(i).pbestNS<gbestNS(it)
                gbest(it,:)=particle(i).pbest;
                gbestNS(it)=particle(i).pbestNS;
            end
        end
    else
        gbest(it,:)=gbest(it-1,:);
        gbestNS(it)=gbestNS(it-1);
        for i=1:npop
            particle(i).velocity=w*particle(i).velocity...
                                +c1*rand*(particle(i).pbest-particle(i).position)...
                                +c2*rand*(gbest(it,:)-particle(i).position);
                            
            particle(i).velocity=min(max(particle(i).velocity,-vmax),vmax);
            
            particle(i).position=particle(i).position+particle(i).velocity;
            
            particle(i).position=min(max(particle(i).position,xmin),xmax);
            
            particle(i).NS=NS(particle(i).position);
            
            if particle(i).NS<particle(i).pbestNS
                particle(i).pbest=particle(i).position;
                particle(i).pbestNS=particle(i).NS;

                if particle(i).pbestNS<gbestNS(it)
                    gbest(it,:)=particle(i).pbest;
                    gbestNS(it)=particle(i).pbestNS;
                end
            end
        end
    end
    
    disp(['Iteration ' num2str(it) ':   Best NS = ' num2str(gbestNS(it))]);
    
    w=w*wdamp;
end
PARSout;
NSout;
plot(gbestNS);

Ntot = sst-1 
ModPred=AFFDEF(gbest(it,:));
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