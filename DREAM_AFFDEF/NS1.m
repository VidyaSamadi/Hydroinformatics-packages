function objfun = NS1(nopt,Pars)
global ss PARSout NSout
load Obs94567.txt
SimRR=AFFDEF(Pars);
objfun=-1*(1-((sum((Obs94567-SimRR).^2))/(sum((Obs94567-mean(Obs94567)).^2))));
NSout(ss,:)=objfun;
PARSout(ss,:)=Pars;
ss=ss+1


