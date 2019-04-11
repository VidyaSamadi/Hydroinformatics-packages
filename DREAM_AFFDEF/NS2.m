function objfunT = NS2(Pars)
global ss PARSout NS2out PopSize Spar Qobs num sst NSt PARt
% Runs the HYMOD model
load Obs94567.txt

Qobs=Obs94567;

for t1=1:PopSize
     pars=Pars(t1,:);
SimRR=AFFDEF(pars);
objfun=-1*(1-((sum((Obs94567-SimRR).^2))/(sum((Obs94567-mean(Obs94567)).^2))));
NS2out(ss,:)=objfun;
PARSout(ss,:)=pars;
ss=ss+1
objfunT(t1,:)=objfun;
if objfun<=-0.5
   Spar(:,sst)=SimRR;
      PARt(sst,:)=pars;
   NSt(sst,:)=objfun;
   sst=sst+1;
end
end
