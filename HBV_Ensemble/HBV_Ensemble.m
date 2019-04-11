%%%------------------------------------------------------------
%%%  Program Name: HBV-Enesemble
%%% 
%%% Citation:
%%% AghaKouchak A., Nakhjiri N., and Habib E., 2013, An educational model for ensemble streamflow 
%%% simulation and uncertainty analysis, Hydrology and Earth System Sciences, 17, 445-452, 
%%% doi:10.5194/hess-17-445-2013.
%%%
%%% AghaKouchak A., Habib E., 2010, Application of a Conceptual Hydrologic Model in Teaching Hydrologic
%%% Processes, International Journal of Engineering Education, 26(4), 963-973.
%%%
%%% Inpt files:
%%% 1) inputPrecipTemp.txt    	(includes 4 columns: (I) date; (II) month ID; (III) Temperature (celcus); and (IV) precipitation (mm/day).
%%% 2) inputMonthlyTempEvap.txt  (includes 3 columns: (I) Monthly Temperature (celcus); (II) Monthly Potential Evapotranspiration; and (III) Mean Daily Potential Evapotranspiration (mm/day).
%%% 3) Qobs.txt            	    (includes 1 column: Observed Runoff (Q) in cubic meters per second (CMS)
%%% 4) BV.txt                    Lower and upper bounds of the model parameters
%%% 5) IV.txt                    Initial conditions, watershed information and number of simulations  
%%%
%%% Note: The parameter estimation and uncertainty analysis is based on the
%%% GLUE approach
%%%
%%% Disclaimer:
%%% This program (hereafter, software) is designed for instructional and educational use only.
%%% Research and commercial use is prohibited. The software is provided 'as is' without warranty 
%%% of any kind, either express or implied. The software could include technical or other mistakes,
%%% inaccuracies or typographical errors. The use of the software is done at your own discretion and 
%%% risk and with agreement that you will be solely responsible for any damage and that the authors
%%% and their affiliate institutions accept no responsibility for errors or omissions in the software
%%% or documentation. In no event shall the authors or their affiliate institutions be liable to you or
%%% any third parties for any special, indirect or consequential damages of any kind, or any damages whatsoever.
%%%
%%% For technical questions contact:
%%% Dr. Amir AghaKouchak
%%% Dept. of Civil & Environmental Engineering
%%% University of California, Irvine
%%% amir.a@uci.edu
%%% 
%%% Acknowledgment:
%%% We are grateful to many colleagues and graduate students who offered valuable comments and suggestions
%%% for improvements. These individuals include Leonardo Valerio Noto, Ali Mehran, Jeff Tuhtan, Nasrin Nasrollahi,
%%% Mehdi Rezaeian Zadeh, Naveen Duggi and Mehdi Javadian. I offer special thanks to students of Watershed
%%% Modeling class (University of California, Irvine) who participated in a survey on their learning gains. 
%%% I apologize to those whose names unintentionally might have been omitted.
%%%
%%%------------------------------------------------------------
clear 
clc

% Reading input file
% read file from inputPrecipTemp.txt
PTIf=fopen('inputPrecipTemp.txt');
% The structure of the file is [date,month1,year,month,Temp,prec]
PTI=textscan(PTIf,'%d/%d/%d %d %f %f'); 
fclose(PTIf);
month_list=PTI{4};
Temp=PTI{5};
prec=PTI{6};


% read file from Qobs.txt
QOf=fopen('Qobs.txt');
qo=textscan(QOf,'%f');
fclose(QOf);
qo=qo{1};

% read file from inputMonthlyTempEvap.txt
MTEf=fopen('inputMonthlyTempEvap.txt');
% The structu
l1u=ub(6); % L ub
k1u=ub(7); % K_1 ub
k2u=ub(8); % K_2 ub
kpu=ub(9); % K_p ub
pwpu=ub(10); % PWP ub

ivf=fopen('IV.txt');
iv=textscan(ivf,'%s %f');
fclose(ivf);

ca=iv{2}(1); % Watershed Area
tt=iv{2}(2); % Snow Melt Thr T_t
ns=iv{2}(3); % number of simulations
snow(1,1:ns)=iv{2}(4); % initial value Snow
soil(1,1:ns)=iv{2}(5); % initial value Soil Moisture 
s1(1,1:ns)=iv{2}(6); % initial value S_1
s2(1,1:ns)=iv{2}(7); % initial value S_2
in = iv{2}(8); % Optimization Criteria 
ensemble=iv{2}(9); % Ensemble mode; "1" is on 

n=size(Temp,1);
aq=mean(qo);
tprec=0;
tea1=0;
reg1=0; %#ok<*NASGU>
nash=0; % Nash Sutcliff 
mqo=sum(qo);

% random generation of values for parameters
d=(rand(ns,1)*(ddu-ddl))+ddl;
fc=(rand(ns,1)*(fcu-fcl))+fcl;
beta=(rand(ns,1)*(betau-betal))+betal;
c=(rand(ns,1)*(cu-cl))+cl;
k0=(rand(ns,1)*(k0u-k0l))+k0l;
l=(rand(ns,1)*(l1u-l1l))+l1l;
k1=(rand(ns,1)*(k1u-k1l))+k1l;
k2=(rand(ns,1)*(k2u-k2l))+k2l;
kp=(rand(ns,1)*(kpu-kpl))+kpl;
pwp=(rand(ns,1)*(pwpu-pwpl))+pwpl;
% end of random generation


% Calculating Q simulated for entered number of simulations
for s=1:ns
    for i=2:n
        if(Temp(i)<tt)
            snow(i,s)=snow(i-1,s)+prec(i);
            lwater(i,s)=0;
            pe(i,s)=(1+c(s)*(Temp(i)-monthly(month_list(i)+1)))*dpem(month_list(i)+1);
            if(soil(i-1,s)>pwp(s))
                ea(i,s)=pe(i,s);
            else
                ea(i,s)=pe(i,s)*(soil(i-1,s)/pwp(s));
            end
            dq(i,s)=lwater(i,s)*((soil(i-1,s)/fc(s))^beta(s));
            soil(i,s)=soil(i-1,s)+lwater(i,s)-dq(i,s)-ea(i,s);
            s1(i,s)=s1(i-1,s)+dq(i,s)-(max(0,s1(i-1,s)-l(s))*k0(s))-(s1(i-1,s)*k1(s))-(s1(i-1,s)*kp(s));
            s2(i,s)=s2(i-1,s)+(s1(i-1,s)*kp(s))-s2(i-1,s)*k2(s);
            q(i,s)=(max(0,s1(i-1,s)-l(s)))*k0(s)+(s1(i,s)*k1(s))+(s2(i,s)*k2(s));
            qm(i,s)=(q(i,s)*ca*1000)/(24*3600);
        else
            snow(i,s)=max(snow(i-1,s)-d(s)*(Temp(i)-tt),0);
            lwater(i,s)=prec(i)+min(snow(i-1,s),d(s)*(Temp(i)-tt));
            pe(i,s)=(1+c(s)*(Temp(i)-monthly(month_list(i)+1)))*dpem(month_list(i)+1);
            if(soil(i-1,s)>pwp(s))
                ea(i,s)=pe(i,s);
            else
                ea(i,s)=pe(i,s)*(soil(i-1,s)/pwp(s));
            end
            dq(i,s)=lwater(i,s)*((soil(i-1,s)/fc(s))^beta(s));
            soil(i,s)=soil(i-1,s)+lwater(i,s)-dq(i,s)-ea(i,s);
            s1(i,s)=s1(i-1,s)+dq(i,s)-(max(0,s1(i-1,s)-l(s))*k0(s))-(s1(i-1,s)*k1(s))-(s1(i-1,s)*kp(s));
            s2(i,s)=s2(i-1,s)+(s1(i-1,s)*kp(s))-s2(i-1,s)*k2(s);
            q(i,s)=(max(0,s1(i-1,s)-l(s)))*k0(s)+(s1(i,s)*k1(s))+(s2(i,s)*k2(s));
            qm(i,s)=(q(i,s)*ca*1000)/(24*3600);
        end
    end
    mq(s,1)=sum(qm(:,s));
    e(:,s)=qo-qm(:,s);% calculating difference between q observed and q measured
    nse(s,1)=1-(sum((e(:,s)).^2))/(sum((qo-aq).^2)); % calculating nse
    qrms(s,1)=sqrt(mean((e(:,s)).^2)); % calculating root mean square value
    regg1(s,1)=corr(qm(:,s),qo);
    tbias(s,1)=mq(s,1)/mqo;
end
nseqm=[nse qm'];

% end of simulation

corr1=regg1;
qrms1=qrms;
nse1=nse;
cb=1-power(tbias,2);
cb1=cb;

% code to get the optimal parameter
if(in==1)
    peako=min(qrms1);%  Calculating minimum of root mean square value to determine optimal parameters
    for p=1:ns
        if(qrms1(p)==peako)
            v=p;% index of optimal parameter
        end
    end
end
if(in==2)
    peako=max(nse1);
    for p=1:ns
        if(nse1(p)==peako)
            v=p;% index of optimal parameter
        end
    end
end
if(in==3)
    peako=max(corr1);
    for p=1:ns
        if(corr1(p)==peako)
            v=p;% index of optimal parameter
        end
    end
end
if(in==4)
    peako=min(cb1);
    for p=1:ns
        if(cb1(p)==peako)
            v=p;% index of optimal parameter
        end
    end
end
cd=d(v);
cfc=fc(v);
cbeta=beta(v);
cc=c(v);
ck0=k0(v);
cl=l(v);
ck1=k1(v);
ck2=k2(v);
ckp=kp(v);
cpwp=pwp(v);
tcq=0;
tqo=0;

csnow(1,1)=iv{2}(4);
csoil(1,1)=iv{2}(5);
cs1(1,1)=iv{2}(6);
cs2(1,1)=iv{2}(7);

% simulation to determine optimal Q
for r=2:n
    tprec=tprec+prec(r);
    if(Temp(r)<tt)
        csnow(r,1)=csnow(r-1,1)+prec(r);
        clwater(r,1)=0;
        cpe(r,1)=(1+cc*(Temp(r)-monthly(month_list(r)+1)))*dpem(month_list(r)+1);
        if(csoil(r-1,1)>cpwp)
            cea(r,1)=cpe(r,1);
        else
            cea(r,1)=cpe(r,1)*(csoil(r-1,1)/cpwp);
        end
        cdq(r,1)=clwater(r)*((csoil(r-1,1)/cfc)^cbeta);
        csoil(r,1)=csoil(r-1,1)+clwater(r,1)-cdq(r,1)-cea(r,1);
        cs1(r,1)=cs1(r-1,1)+cdq(r,1)-(max(0,cs1(r-1,1)-cl)*ck0)-(cs1(r-1,1)*ck1)-(cs1(r-1,1)*ckp);
        cs2(r,1)=cs2(r-1,1)+(cs1(r-1,1)*ckp)-cs2(r-1,1)*ck2;
        cq(r,1)=(max(0,cs1(r-1,1)-cl))*ck0+(cs1(r,1)*ck1)+(cs2(r,1)*ck2);
        cqm(r,1)=(cq(r,1)*ca*1000)/(24*3600);
        cqs(r,1)=((qo(r,1)-cqm(r,1))^2);
        cqms(r,1)=((qo(r,1)-aq)^2);
    else
        csnow(r,1)=max(csnow(r-1,1)-cd*(Temp(r)-tt),0);
        clwater(r,1)=prec(r)+min(csnow(r-1,1),cd*(Temp(r)-tt));
        cpe(r,1)=(1+cc*(Temp(r)-monthly(month_list(r)+1)))*dpem(month_list(r)+1);
        if(csoil(r-1,1)>cpwp)
            cea(r,1)=cpe(r,1);
        else
            cea(r,1)=cpe(r,1)*(csoil(r-1,1)/cpwp);
        end
        cdq(r,1)=clwater(r)*((csoil(r-1,1)/cfc)^cbeta);
        csoil(r,1)=csoil(r-1,1)+clwater(r,1)-cdq(r,1)-cea(r,1);
        cs1(r,1)=cs1(r-1,1)+cdq(r,1)-(max(0,cs1(r-1,1)-cl)*ck0)-(cs1(r-1,1)*ck1)-(cs1(r-1,1)*ckp);
        cs2(r,1)=cs2(r-1,1)+(cs1(r-1,1)*ckp)-cs2(r-1,1)*ck2;
        cq(r,1)=(max(0,cs1(r-1,1)-cl))*ck0+(cs1(r,1)*ck1)+(cs2(r,1)*ck2);
        cqm(r,1)=(cq(r,1)*ca*1000)/(24*3600);
        cqs(r,1)=((qo(r,1)-cqm(r,1))^2);
        cqms(r,1)=((qo(r,1)-aq)^2);
    end
    tea1=tea1+cea(r,1);
    tcq=tcq+cqm(r,1);
end
% end of simulation

e1=qo-cqm;% calculating difference between q observed and q measured
cnse=1-(sum(e1.^2))/(sum((qo-aq).^2));
tdis=tprec-tea1;
reg1=corr(cqm,qo);
bias=((tcq-mqo)/mqo)*100;

%Set the results
result.snow=csnow;
result.lwater=clwater;
result.soil=csoil;
result.dq=cdq;
result.s1=cs1;
result.pe=cpe;
result.ea=cea;
result.s2=cs2;
result.qsim=cq;
result.correlation=reg1;
result.nash=cnse;
result.mp=[cd;cfc;cbeta;cc;ck0;cl;ck1;ck2;ckp;cpwp];


figure;
plot(cqm,'k-','linewidth',1);
hold on
plot(qo,'r','linewidth',1);
hold off
h = legend('Simulated','Observed',1);
set(h,'Interpreter','none')
xlabel('Time');
ylabel('Runoff');
ylmam1=max(cqm);
ylmax2=max(qo);
ymaxs=[ylmam1 ylmax2];
ylim([0 max(ymaxs)])

Qsim=cqm;
save Qsim.mat Qsim

% For Ensemble mode------
if ensemble==1
    qmAccept=[];

    for i=1:size(nseqm,1)
        if nseqm(i,1)>0.5
            qmAccept=[qmAccept; nseqm(i,2:size(nseqm,2))];
        end
    end
    QsimALL=qmAccept;

    figure;
    plot(cqm,'k-','linewidth',1);
    hold on
    plot(qo,'r','linewidth',1);
    for jk=1:size(qmAccept,1)
        plot(qmAccept(jk,2:size(qmAccept,2)),'color',[0.7,0.7,0.7],'linewidth',2);
        hold on
    end
    plot(cqm,'k-','linewidth',1);
    hold on
    plot(qo,'r','linewidth',1);
    hold off
    h = legend('Simulated','Observed','Uncertainity',1);
    set(h,'Interpreter','none')
    titlestr = 'Ensemble Criteria: Nash-Sutcliffe Coefficient Above 0.5';
    title(titlestr)
    
    save QsimALL.mat QsimALL
end