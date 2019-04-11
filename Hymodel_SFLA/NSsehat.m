
% Runs the HYMOD model
load bound.txt;
num=1850;
Qobs=bound(1:num,4);
PET=bound(1:num,5);
Precip=sum(bound(1:num,6:9),2);
% Define the rainfall
MaxT =num;
% Define the parameters
cmax = Pars(1); bexp = Pars(2); alpha = Pars(3); Rs = Pars(4); Rq = Pars(5);
% HYMOD PROGRAM IS SIMPLE RAINFALL RUNOFF MODEL
x_lobound = 0.0;
% Initialize slow tank state
x_slow = 2.3503/(Rs*22.5);
% Initialize state(s) of quick tank(s)
x_quick(1:3,1) = 0; t=1; outflow = [];
% START PROGRAMMING LOOP WITH DETERMINING RAINFALL - RUNOFF AMOUNTS
while t<MaxT+1,
   Pval = Precip(t,1); PETval = PET(t,1);
   % Compute excebound precipitation and evaporation
   [UT1,UT2,x_lobound] = excess(x_lobound,cmax,bexp,Pval,PETval);
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
objfun=(1-((sum((bound(1:num,4)-SimRR).^2))/(sum((bound(1:num,4)-mean(bound(1:num,4))).^2))))
