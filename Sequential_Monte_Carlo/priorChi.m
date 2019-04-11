%% Prior : Inverse Chi distribution
function[varargout]=priorChi(v,l,tita)

prior=((v/2)*log((l*v)/2))-gammaln(v/2)+((-v*l)/(2*tita))-((1+(v/2))*log(tita));

varargout{1}=prior;




