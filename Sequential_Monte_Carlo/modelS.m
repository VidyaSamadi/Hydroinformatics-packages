%% Identifying the initial mean and variance for the C parameters 
function[varargout]=modelS(f,e)

wei_mean=e*(gamma(1+(1/f)));
wei_var =(e^2*(gamma(1+(2/f))))-(wei_mean)^2;

varargout{1}=wei_mean;
varargout{2}=wei_var;