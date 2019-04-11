%% Prior : Weibull distribution
function[varargout]=priorWeibull(a,lam,tita)

prior=(a/lam)*((tita/lam)^(a-1))*exp(-(tita/lam)^a);

varargout{1}=prior;




