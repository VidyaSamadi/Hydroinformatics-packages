%% Identifying the initial mean and variance for the parameter variance 
function[varargout]=modelV(v,l)

chi_mean=(v*l)/(v-2);
chi_var=(2*v^2*l^2)/(((v-2)^2)*(v-4));

varargout{1}=chi_mean;
varargout{2}=chi_var;