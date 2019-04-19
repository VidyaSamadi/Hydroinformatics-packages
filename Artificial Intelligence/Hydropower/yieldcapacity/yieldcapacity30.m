clear
clc
Idata = xlsread('data');
I=reshape(Idata,360,1);
EnvDemand=10000;
Smin=2000;
Alphatol=0.01;
Deltay=20;
for Alphastr=[0.75 0.85 0.9]
    
for j=5:55;    
    K(j)=6000*j;
    
    Sini(j)=K(j);
         
    R=zeros(360,1);
    Z=zeros(360,1);
    S=zeros(360,1);
    i=0;
    stop=0;
    
    while stop==0
        i=i+1;
        if i==1
            y(i)=min(I);
        end
        S(1)=Sini(j);
        
        for t=1:360;
            if S(t)+I(t)<EnvDemand+Smin
                R(t)=0;
                Z(t)=0;
                S(t+1)=Smin;
            elseif S(t)+I(t)>=EnvDemand+y(i)+Smin && S(t)+I(t)<EnvDemand+y(i)+K(j)
                R(t)= y(i);
                Z(t)=1;
                S(t+1)=S(t)+I(t)-R(t)-EnvDemand;
            elseif S(t)+I(t)>= EnvDemand+y(i)+K(j)
                R(t)=I(t)+S(t)-K(j)-EnvDemand ;
                S(t+1)=K(j);
                Z(t)=1;
            else
                R(t)= S(t)+I(t)-EnvDemand-Smin;
                Z(t)=0;
                S(t+1)=Smin;
            end
            
            
        end
        alpha(i)=mean(Z);
        if alpha(i)< Alphastr - Alphatol
            
            y(i+1)=y(i)-Deltay;
        elseif alpha(i)> Alphastr+Alphatol
            y(i+1)=y(i)+Deltay;
        else
            stop=1;
            ystr(j)=y(i);
            if ystr(j)<0
            ystr(j)=0;
            end
        end
    end
end
    plot(K,ystr);
    hold on;
    
    
end

ylabel('Yield')
xlabel('Ka')
title('Safe Yield')
