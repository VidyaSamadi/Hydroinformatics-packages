function [p] = normalfunc(x,Extra)
% Normal pdf

% MATLAB function
%p = mvnpdf(x,Extra.mu,Extra.qcov);

% Log density
p = -0.5*x*Extra.invC*x';