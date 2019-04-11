
load Obs94567.txt
SimRR=AFFDEF(Pars);
NS=(1-((sum((Obs94567-SimRR).^2))/(sum((Obs94567-mean(Obs94567)).^2))));

