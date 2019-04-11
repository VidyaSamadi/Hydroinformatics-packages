
load bound.txt
global num
num=1850;
Qobs=bound(1:num,4);
ModPred=SimRR(Px);
% Compute the RMSE
RMSE = sqrt ( sum ( ( ModPred - Qobs).^2) / num);
for i=1:Ntot
Spar(:,i)=SimRR(pars(i,:));
end

% Now create the parameter and total uncertainty
for j = 1:Ntot
    Smod(:,j) = Spar(:,j) + normrnd(0,RMSE,num,1);
end;



% Now sort to get 95% ranges
for j = 1:num
    a = sort(Spar(j,1:Ntot));
    % And take the 2.5 and 97.5 percentiles
    minpar95(j,1) = a(floor(0.025 * Ntot)); maxpar95(j,1) = a(floor(0.975 * Ntot));
    % Same with total uncertainty
    a = sort(Smod(j,1:Ntot));
    % And take the 2.5 and 97.5 percentiles
    mintot95(j,1) = a(floor(0.025 * Ntot)); maxtot95(j,1) = a(floor(0.975 * Ntot));
end;

% Now plot the results 
plot([1:num],minpar95,[1:num],maxpar95); % These are the upper and lower parameter uncertainty bounds
hold on; plot([1:num],mintot95,'R',[1:num],maxtot95,'G'); % These are the upper and lower total uncertainty bounds
plot([1:num],Qobs,'r+');

% Now calculate percentage inside the 95% bound --> should be about 95%
Contained = 1 - length(find(Qobs<mintot95 |  Qobs> maxtot95))/num
Contained1 = 1 - length(find(Qobs<minpar95 |  Qobs> maxpar95))/num
toc
% Note the ranges do not look very good --> because of homoscedasticity. If
% the error is assumed to increase with flow level results look much better

