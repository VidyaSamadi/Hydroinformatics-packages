
clc;
clear;

npop=5;
nvar=1;

w=1;
wdamp=0.99;

c1=2;
c2=2;

xmin=6000;
xmax=100000000;
dx=xmax-xmin;

vmax=0.1*dx;

maxit=100;

empty_particle.position=[];
empty_particle.velocity=[];
empty_particle.Cost2=[];
empty_particle.pbest=[];
empty_particle.pbestCost2=[];

particle=repmat(empty_particle,npop,1);

gbest=zeros(maxit,nvar);
gbestCost2=zeros(maxit,1);

for it=1:maxit
    if it==1
        gbestCost2(1)=inf;
        for i=1:npop
            particle(i).velocity=zeros(1,nvar);
            particle(i).position=xmin+(xmax-xmin)*rand(1,nvar);
            particle(i).Cost2=Cost2(particle(i).position);
            
            particle(i).pbest=particle(i).position;
            particle(i).pbestCost2=particle(i).Cost2;
            
            if particle(i).pbestCost2<gbestCost2(it)
                gbest(it,:)=particle(i).pbest;
                gbestCost2(it)=particle(i).pbestCost2;
            end
        end
    else
        gbest(it,:)=gbest(it-1,:);
        gbestCost2(it)=gbestCost2(it-1);
        for i=1:npop
            particle(i).velocity=w*particle(i).velocity...
                                +c1*rand*(particle(i).pbest-particle(i).position)...
                                +c2*rand*(gbest(it,:)-particle(i).position);
                            
            particle(i).velocity=min(max(particle(i).velocity,-vmax),vmax);
            
            particle(i).position=particle(i).position+particle(i).velocity;
            
            particle(i).position=min(max(particle(i).position,xmin),xmax);
            
            particle(i).Cost2=Cost2(particle(i).position);
            
            if particle(i).Cost2<particle(i).pbestCost2
                particle(i).pbest=particle(i).position;
                particle(i).pbestCost2=particle(i).Cost2;

                if particle(i).pbestCost2<gbestCost2(it)
                    gbest(it,:)=particle(i).pbest;
                    gbestCost2(it)=particle(i).pbestCost2;
                end
            end
        end
    end
    
    disp(['Iteration ' num2str(it) ':   Best Cost2 = ' num2str(gbestCost2(it))]);
    
    w=w*wdamp;
	plot(gbestCost2(1:it));

	pause(0.1)
end

plot(gbestCost2);
title('convergence')

disp(['Best Ka=',num2str(gbest(end)')]);


