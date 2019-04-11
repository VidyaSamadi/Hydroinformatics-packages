function z=fitness(pars)
    num=2557;
    load Hisflow25.txt
    Hisflow25=Hisflow25(1:2557);
    SimRR=ARNO(pars);
   % ND = find(Hisflow25(1:num)>=4);
    % Find driven part of hydrograph
   % D = find(Hisflow25(1:num)<1);
    
ss=0;
ff=0;

for i=1:num
 ss=ss+((1/2557)*((log10(Hisflow25(i)/ SimRR(i))).^2));
end
for j=1:num
  ff=ff+((1/2557)*(sum((Hisflow25(j)-SimRR(j)).^2)));
end
 
    

 z1=sqrt(ss);
    z2=sqrt(ff);
 %z3=((1/2557)*(abs(sum(Hisflow25) - sum(SimRR))/sum(Hisflow25)));
 
    z=[z1 z2]';

end