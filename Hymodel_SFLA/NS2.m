function objfunT = NS2(Pars)
global ss PARSout NS2out PopSize Spar Qobs num sst NSt PARt
% Runs the HYMOD model
load bound.txt;

Qobs=bound(1:num,4);
PET=bound(1:num,5);
Precip=sum(bound(1:num,6:9),2);
% Define the rainfall
MaxT =num;
for t1=1:PopSize
     pars=Pars(t1,:);
% Define the parameters
cmax = pars(1); bexp = pars(2); alpha = pars(3); Rs = pars(4); Rq = pars(5);
% HYMOD PROGRAM IS SIMPLE RAINFALL RUNOFF MODEL
x_loss = 0.0;
% Initialize slow tank state
x_slow = 2.3503/(Rs*22.5);
% Initialize state(s) of quick tank(s)
x_quick(1:3,1) = 0; t=1; outflow = [];
% START PROGRAMMING LOOP WITH DETERMINING RAINFALL - RUNOFF AMOUNTS
while t<MaxT+1,
   Pval = Precip(t,1); PETval = PET(t,1);
   % Compute excess precipitation and evaporation
   [UT1,UT2,x_loss] = excess(x_loss,cmax,bexp,Pval,PETval);
   % Partition UT1 and UT2 into quick and slow flow component
   UQ = alpha*UT2 + UT1; US = (1-alpha)*UT2;
   % Route slow flow component with single linear reservoir
   inflow = US; [x_slow,outflow] = linres(x_slow,inflow,outflow,Rs); QS = outflow;
   % Route quick flow component with linear reservoirs
   inflow = UQ; k = 1; 
   while k < 4,
      [x_quick(k),outflow] = linres(x_quick(k),inflow,outflow,Rq); inflow = outflow; 
      k = k+1;
   end;
   % Compute total flow for timestep
   output(t,1) = (QS + outflow)*22.5;
   t = t+1;   

end;
   SimRR = output(1:MaxT,1);
objfun=-1*(1-((sum((bound(1:num,4)-SimRR).^2))/(sum((bound(1:num,4)-mean(bound(1:num,4))).^2))));
NS2out(ss,:)=objfun;
PARSout(ss,:)=pars;
ss=ss+1
objfunT(t1,:)=objfun;
if objfun<=-0.65
   Spar(:,sst)=SimRR;
      PARt(sst,:)=pars;
   NSt(sst,:)=objfun;
   sst=sst+1;
end
end
