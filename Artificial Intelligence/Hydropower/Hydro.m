clear 
close all
clc

% Assign Initial Values
load IData.mat              % Read Inflow Data (in CMS)

IData=IData*30*24*3600/(10^6); % Convert Discharge(m3/S) to Volume(MCM)

% Reservoir Inflow Matrix
[m,n]=size(IData);

I=reshape(IData',m*n,1);% Convert Matrix to Vector

T=m*n;

syms x;

Vol=0.0008*(x^3)-1.8452*(x^2)+1348.2*x+328603 ; %Elevation-Volume Curve

hd=52 ;                  % Design Head(m)																			!

htail=722;% Tail Water Elevation(m)

hloss=2; %Loss Head(m)

hmax=1.25*hd ;%Maximum Elevation(m)


hmin=0.65*hd ;%Minimum Elevation(m)
Qd=435.6*30*24*3600*(10^-6); %Design Discharge(MCM)
NU=1; % Number of Unit
Pmax=50*NU; % Plant Capacity(MW)
Emax=Pmax*30*24;% Plant Capacity(MWh)
Qmax=1.1*Qd ;%Maximum Discharge
Qmin=0.5*Qd/NU ;%Minimum Discharge
e=0.8 ; %Generation Efficiency
AlphaTol=0.01; %Delta Alpha
DeltaFE=100; %Delta Firm Energy
MOL=279;%Minimum Operation Level(m)
MOV=double(subs(Vol,MOL));%Minimum Operation Volume(MCM)
ss=MOV ;%Initial Volume

AlphaStar= 0.9; % Target reliability
    display(['Run for AlphaStar=' num2str(AlphaStar),'  >>>>>Please Waite...'])
for ka=[1 200 500 800 1000 1200 1600 1800 2000 2200];% Active Storage
if (ka+MOV)>1707 %1707 is Storage(MCM) that It's Level is Equal to Maximum Height for Energy Generation  
    M=1707;
elseif (ka+MOV)<=1707 
    M=ka+MOV;
end
%% Estimation of Firm Energy (FE)
%Pre allocation
hh=zeros(1,T);
RR=zeros(1,T);
hhbar=zeros(1,T);
hhnet=zeros(1,T);
ssbar=zeros(1,T);
HE=zeros(1,T);
RE=zeros(1,T);
E=zeros(1,T);
z=zeros(1,T);
j=0; % Iteration Index
stop=0;
while stop==0
    j=j+1;
    if j==1
        FE(j)=2.72*min(I)*((hmax+hmin)/2)*e; % First Guess of FE
    end % end of if
    for t=1:T
        % Estimation of R(t) for FE(j)
        k=1;
        s=[50 100];
       while  abs(s(k)-s(k+1))>10^-5 
             k=k+1; 
            if k==2
                s(k)=ss(t);
            end
                sbar(k)=(ss(t)+s(k))/2;%Mean Storage
                v=[0.7292 -375.92 +48424-sbar(k)];
                hbar(k)=max(roots(v));
                hnet(k)=hbar(k)-htail-hloss;% Net Head
                R(k)=FE(j)/(2.72*hnet(k)*e);% Release
                s(k+1)=ss(t)+I(t)-R(k);    
       end % End of Second while
           ss(t+1)=s(k+1);
           RR(t)=R(k);      
       if  s(k+1)<MOV
           ss(t+1)=MOV;
           RR(t)=R(k)-(MOV-s(k+1));
       elseif s(k+1)>M
           ss(t+1)=M;
           RR(t)=R(k)+(s(k+1)-M);
       end %End of if
       ssbar(t)=(ss(t)+ss(t+1))/2;%Mean Storage
       vt=[0.7292 -375.92 +48424-ssbar(t)];
       hhbar(t)=max(roots(vt));
       hhnet(t)=hhbar(t)-htail-hloss;%Net Head 
       HE(t)=hhnet(t);
       RE(t)=RR(t);
       if HE(t)<hmin || HE(t)>hmax
            HE(t)=0;
       end
%        if RE(t)<Qmin 
%             RE(t)=0;
%        end%End of if
%        if RE(t)>Qmax 
%             RE(t)=Qmax;
%        end
        E(t)=2.72*HE(t)*RE(t)*e;% Energy Calculation
        if E(t)>Emax 
            E(t)=Emax;
        end%End of if
        if E(t)>=FE(j)
            z(t)=1;
        else
            z(t)=0;
        end % End of if      
    end%End of if      
Alpha(j)=sum(z)/T;
if Alpha(j)<(AlphaStar-AlphaTol)
    FE(j+1)=FE(j)-DeltaFE;
elseif Alpha(j)>(AlphaStar+AlphaTol)
    FE(j+1)=FE(j)+DeltaFE;
else
    FEStar=FE(j);
    stop=1;
end %End of if
end %End of First while
%% Save Output
reslt(ka,(AlphaStar*100))=ka;
reslt(ka,(AlphaStar*100)+1)=FE(j);
reslt(ka,(AlphaStar*100)+2)=AlphaStar;
end %End of for ka


for i=numel(reslt(:,1)):-1:1
    if reslt(i,:)==0
        reslt(i,:)=[];
    end %End of if
end %End of for
for i=numel(reslt(1,:)):-1:1
    if reslt(:,i)==0
        reslt(:,i)=[];
    end %End of if
end %End of for
display('END of RUN')
