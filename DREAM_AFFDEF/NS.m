function objfun = NS(nopt,Pars)
global ss PARSout NSout Spar Qobs num sst NSt PARt
% Runs the HYMOD model
load Obs94567.txt
Qobs=Obs94567;

SimRR=AFFDEF(Pars);
objfun=-1*(1-((sum((Obs94567-SimRR).^2))/(sum((Obs94567-mean(Obs94567)).^2))));
NSout(ss,:)=objfun;
PARSout(ss,:)=Pars;
ss=ss+1
if objfun<=-0.5
   Spar(:,sst)=SimRR;
      PARt(sst,:)=Pars;
   NSt(sst,:)=objfun;
   sst=sst+1;
end

