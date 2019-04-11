%% Identifying the initial mean and variance for the K and BFI parameters 
function[varargout]=modelMB(tao,a,b,x,y)

beta_mean=a/(a+b);
beta_var=(a*b)/((a+b)^2*(a+b+1));
beta_mean2=x/(x+y);
beta_var2=(x*y)/((x+y)^2*(x+y+1));
mbeta_mean=(tao*beta_mean)+((1-tao)*(beta_mean2));
mbeta_var = (tao*beta_var)+((1-tao)*(beta_var2));

varargout{1}=mbeta_mean;
varargout{2}=mbeta_var;