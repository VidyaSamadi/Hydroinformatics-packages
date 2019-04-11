function f=Hedging(Kp,Extra)
%Hedging code with Kp=1 (SOP)
Kp=[12.1423356549441,11.5919553349398,11.2169924375210,9.34410403629298,7.81078738313953,8.96228261341838,6.42261465827490,1.22966163235605,1.05365314566976,12.3930132692607,1.08844275759474,11.2539672804274];
%St1=Kp(13);
%Kp=Kp(1:12);
%Kp=ones(1,12); 
%Kp=[1.05154391918954,1.75951597953507,1.27537229747578,1.04740370110532,18.3405137795838,20,2.56537091121742,1.62565075608919,1,1,1.30541537924786,1];
global NK z2 ss st100 Send kp RT UTprc ST NumRun
NumRun=NumRun+1
z2=0;
NK=[0 0]; %number of months with NK(1)= is related to Necessary Demand (Recharge,2.66, and drink, Sh(i), and Evaporation E(i)) and NK(2) concerned secondaty Demand, Plain Qazvin
Send=[];
C=320; %Reservoir active Capacity

%it=[33.2	61	34.3	9.9	4.5	3.5	3.1	3.3	3.5	3.5	4.4	14.1
%22.6	62.3	50.8	15.5	6.2	4.2	8.5	10.3	5.6	4.8	6.5	13.2
%22.3	60	39.4	14.8	6.8	4	6.4	6.2	4.3	4.5	4.1	6.1
%12	29.8	33.4	9.2	3.6	3.4	3.2	3.4	3.6	3.2	3.6	10.1
%22.7	39.2	53.7	25.7	7.4	4.8	4.1	11.6	9.9	8.4	7.3	29.2
%95.9	154.9	101.2	48.4	19	10.7	9.8	15.5	7.1	6.5	6.3	8.1
%29.5	26.5	25	9.5	5.5	3.7	3.6	3.5	3.9	4.4	4.1	10
%26.1	64.3	46.5	16.1	6.4	4	7.2	5.4	6.4	4.3	4.3	8.4
%32.2	61.3	71.9	34.5	11.2	6.1	5.4	8.8	6.2	5.7	6.5	14.7
%23.8	50.1	41.8	15.2	6.8	5.3	3.6	4.4	4.4	3.6	3.6	9.8
%29.1	35.5	30.4	13.1	4.6	6.5	4.3	3.7	3.5	3.7	3.2	7.4
%33.2	50.4	42.6	16.3	6.4	4.4	4.1	3.8	4.3	4.5	4.8	5.5
%28.2	44.7	43.4	22.3	8.6	4.8	4	4.6	4.6	4.4	3.9	9.5
%14.5	23.1	31.8	9.4	4.3	2.6	3.9	7.2	5.9	5.6	5.6	11
%27.6	34.9	28.3	12.6	5.5	4.3	3.8	4.6	6.3	5.4	6.3	8.5
%41.5	46.7	45.3	19.9	8.7	4.8	4.5	6	4.2	3.9	3.9	6
%42.1	38	27.3	12.1	6.1	3.3	4	8	5.8	5.2	6.4	12.5
%30.5	60.2	48.7	30	9.4	5	6	5.3	4.5	3.9	4.2	6.6
%30.1	43.3	21.7	10.1	5.1	3.3	8.8	7.6	6.3	5	4.9	7.2
%19.1	46.7	33.9	12.6	6	4.6	4	4.6	6.7	4.6	5.6	9.2
%22.8	36.6	37.3	15	6.5	4.5	6	7	10.1	6.5	9.1	10.6
%36.3	46.6	36.3	17.3	7.8	6.4	4	4.5	6	5.5	5.6	10.6
%20.6	48.1	41.7	16.7	7.4	4.2	4.4	5.4	5	5.4	7.4	12.8
%35.5	69.9	47.9	20.3	10.1	4.7	6.6	14.8	8.9	9.5	8.3	14.7
%42.6	75.1	46.7	25.1	11	7.1	6.2	6.3	5.9	5.2	5.4	10.7
%26.4	29.6	20.7	9.2	5.1	3.2	3.3	3.6	4	4.1	4.5	9.4
%18.9	36.4	21.5	9.6	4.7	3.3	3.4	5.6	3.4	3.3	3.9	6.5
%24.4	29.5	19.1	8.6	3.9	2.9	3.5	3.9	4.3	4.1	3.3	6.4
%28.1	58.1	72.6	38.5	13	6.2	5.6	5.6	4.9	4.5	4.7	7.2
%24.2	33.4	35.5	15.4	10.3	4.2	4.4	9.1	10.1	11.4	7.8	12
%37.5	62.7	44.3	21.4	9.6	7.3	6	18.9	25.2	10.7	10.1	13.9
%29.8	50.9	43.9	21.8	9.6	5.4	5.2	5.7	4	3.4	5	8.4
%30.9	62.9	39.4	16.9	5.8	4.5	4.5	4.3	3.3	3.4	3.6	5.1
%15.4	43.8	29.9	12.9	4.1	3.3	3.6	5.2	4.5	4.6	6	11.5
%44.6	62.5	40.6	17.8	7.9	4.8	4.8	4.7	4.4	4.3	4.2	6.8
%];
    for i=1:12
        %interpolation for calculating of evaporation demand ,E (MCM)
        if It(j,i)+St((j-1)*12+i)+100>420
            E(i)=A(12)*ET(i)/1000; %calculation of Evaporation loss (E) at time i
        else if It(j,i)+St((j-1)*12+i)+100<volume(1);
                area(i)=A(1);
                E(i)=area(i)*ET(i)/1000;
            else
                area(i)=interp1(volume,A,(It(j,i)+100+St((j-1)*12+i)));
                E(i)=area(i)*ET(i)/1000;
            end
        end
        
        WA(i)=It(j,i)+St((j-1)*12+i)-E(i); % Water Availability at time i 
        
        if WA(i)>=Kp(i)*DT(i)
            if WA(i)-(DT(i)+C)>=0
                Rt(i)=DT(i);
                Yt(i)=1/Kp(i)*(It(j,i)+St((j-1)*12+i))-Rt(i);
                St((j-1)*12+i+1)=C;
            else
                Rt(i)=DT(i);
                Yt(i)=1/Kp(i)*(It(j,i)+St((j-1)*12+i))-Rt(i);
                St((j-1)*12+i+1)=WA(i)-Rt(i);
            end
            Ut(1:2,i)=0;
        end

        
        if  WA(i)<Kp(i)*DT(i)
            Rt(i)=1/Kp(i)*WA(i);
            Yt(i)=0;
            St((j-1)*12+i+1)=WA(i)-Rt(i);
            if  St((j-1)*12+i+1)>C
                %Rt(i)=DT(i);
                St((j-1)*12+i+1)=C;
                W(i)=WA(i)-Rt(i)-C; %Spill in period i
            end
            if Rt(i)>=DT(i)    %necessary demand is supported. Secondary Demand is NOT...
                Ut(:,i)=0;
            else if Rt(i)>Sh(i)+2.66
                    NK(2)=NK(2)+1;
                    Ut(1,i)=0; %Deficit for necessary demand
                    Ut(2,i)=DT(i)-Rt(i); %Deficit for secondary demand
                else  %  necessary demand is not supported. Secondary Demand is NOT...
                    NK(1)=NK(1)+1;
                    NK(2)=NK(2)+1;
                    Ut(1,i)=Sh(i)+2.66-Rt(i);
                    Ut(2,i)=Dt(i);
                end
            end
        end
        
    end
    for p=1:12
        Utprc(p)=sum(Ut(:,p))/DT(p);
    end
    Def(j,:)=sum(Ut,1);
    Ytt(j,:)=Yt; %Slack variable for three years
    M(j,:)=Utprc; 
    RTT(j,:)=Rt;
    UTprc(j,1:12)=Utprc;
    STT(j,:)=St(1,(j-1)*12+1:j*12); % Storage matrix for three years
end
diff=St(1)-St((j-1)*12+i+1);
Send(ss,:)=diff;
kp(ss,:)=Kp(1:12);

f1=sum(sum(M.^2))+abs(diff/10);
f2=1-(NK(2)/(3*12)); %calculation of reliability, minus(-) was applied for optimization setting


f=[f1 f2];

%z=max(M,[],1);
%f=z2+max(z,[],2);
%diff=(abs(St(1)-St((j-1)*12+i+1)))/200;
%f=z2%+diff;
end