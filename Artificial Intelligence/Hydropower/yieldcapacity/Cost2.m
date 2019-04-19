function [C]=Cost2(K)

load Idata.mat
I=reshape(Idata,360,1);
EnvDemand=10000;
Smin=2000;
UrbDemand=15000;
Deltay=20;
 Penalty=0;   
   

    
    Sini=K;
         
    R=zeros(360,1);
    Z=zeros(360,1);
    S=zeros(360,1);
    
    S(1)=Sini;
        
        for t=1:360;   %645
            W(t)=S(t)+I(t)-20000;
            
            if W(t)<0
                Penalty(t)=1;
            end
             
            
            if W(t)<=0.2*EnvDemand+Smin
                R(t)=0.8*UrbDemand;
                EnvFlow(t)=W(t)-R(t)-Smin;
                S(t+1)=Smin;
                Loss(t)= 3000*(EnvDemand-EnvFlow(t))+ 1800*0.8*UrbDemand;
            elseif W(t)> 0.2*EnvDemand+0.2*UrbDemand +Smin && W(t)<0.2*EnvDemand+0.2*UrbDemand+K
                R(t)= UrbDemand;
                EnvFlow(t)=EnvDemand;
                S(t+1)=W(t)-R(t)-EnvDemand;
                Loss(t)= 0;
                Benefit(t)=645*R(t);
            elseif W(t)>= 0.2*EnvDemand+0.2*UrbDemand+K
                EnvFlow(t)=EnvDemand;
                R(t)=W(t)-K-EnvDemand ;
                S(t+1)=K;
                Benefit(t)=645*UrbDemand;
                Loss(t)=0;
            else
                EnvFlow(t)=EnvDemand;
                R(t)= S(t)+I(t)-EnvDemand-Smin;
                S(t+1)=Smin;
                Loss(t)= 1800*(UrbDemand-R(t));
                Benefit(t)=645*R(t);
            end
            
        end
        C=sum(Loss)+2000*K-sum(Benefit)+sum(Penalty)*10000000000;
        
    
end




