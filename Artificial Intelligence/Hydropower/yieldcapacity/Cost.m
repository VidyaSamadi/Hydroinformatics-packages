function [C]=Cost(K)

load Idata.mat
I=reshape(Idata,360,1);
EnvDemand=10000;
Smin=2000;
UrbDemand=15000;
Deltay=20;
    
    
    
    
    Sini=K;
         
    R=zeros(360,1);
    Z=zeros(360,1);
    S=zeros(360,1);
    
    S(1)=Sini;
        
        for t=1:360;   %645
            if S(t)+I(t)<=EnvDemand+Smin
                R(t)=0;
                EnvFlow(t)=S(t)+I(t)-Smin;
                S(t+1)=Smin;
                Loss(t)= 3000*(EnvDemand-EnvFlow(t))+ 1800*UrbDemand;
            elseif S(t)+I(t)>=EnvDemand+UrbDemand +Smin && S(t)+I(t)<EnvDemand+UrbDemand+K
                R(t)= UrbDemand;
                EnvFlow(t)=EnvDemand;
                S(t+1)=S(t)+I(t)-R(t)-EnvDemand;
                Loss(t)= 0;
                Benefit(t)=645*R(t);
            elseif S(t)+I(t)>= EnvDemand+UrbDemand+K
                EnvFlow(t)=EnvDemand;
                R(t)=S(t)+I(t)-K-EnvDemand ;
                S(t+1)=K;
                Benefit(t)=UrbDemand;
                Loss(t)=0;
            else
                EnvFlow(t)=EnvDemand;
                R(t)= S(t)+I(t)-EnvDemand-Smin;
                S(t+1)=Smin;
                Loss(t)= 1800*(UrbDemand-R(t));
                Benefit(t)=645*R(t);
            end
            
        end
        C=sum(Loss)+2000*K-sum(Benefit);
        
    




