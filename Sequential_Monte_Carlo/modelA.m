%% Identifying the initial mean and variance for the A parameters 
function[varargout]=modelA(a,b)

beta_mean=a/(a+b);
beta_var=(a*b)/((a+b)^2*(a+b+1));

varargout{1}=beta_mean;
varargout{2}=beta_var;

