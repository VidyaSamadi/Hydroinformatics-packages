%% Prior : Beta distribution
function[varargout]=priorBeta(a,b,tita)

prior=(gamma(a+b)/(gamma(a)*gamma(b)))*(tita^(a-1))*((1-tita)^(b-1));

varargout{1}=prior;


