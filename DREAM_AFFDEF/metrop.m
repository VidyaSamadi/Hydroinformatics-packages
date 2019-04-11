function [newgen,alpha,accept] = metrop(x,p_x,log_p_x,x_old,p_old,log_p_old,N,Sigma,MCMCPar,Extra,option);
% Metropolis rule for acceptance or rejection

% Calculate the number of Chains
NrChains = size(x,1);

% First set newgen to the old positions in X
newgen = [x_old p_old log_p_old];

% And initialize accept with zeros
accept = zeros(NrChains,1);

% -------------------- Now check the various options ----------------------
if option == 1,
    alpha = min(p_x./p_old,1);
end;

if option == 2 | option == 4 % Lnp probability evaluation
    alpha = min(exp(p_x - p_old),1);
end;

if option == 3, % SSR probability evaluation
    alpha = min((p_x./p_old).^(-N.*(1+MCMCPar.Gamma)./2),1);
end;

if option == 5, % SSR probability evaluation, but now weighted with mesurement error
    % Note that measurement error is single number --> homoscedastic; variance can be taken out of sum sign
    alpha = min(exp(-0.5*(-p_x + p_old)./Sigma^2),1);   % signs are different because we write -SSR
end;

if option == 6, % SSR probability evaluation, but now weighted with mesurement error
    % Note that measurement error is a vector --> heteroscedastic; variance within sum sign  -- see CompDensity.m
    alpha = min(exp(-0.5*(-p_x + p_old)),1);            % signs are different because we write -SSR
end;
% -------------------------------------------------------------------------

% Generate random numbers
Z = rand(NrChains,1);

% Find which alpha's are greater than Z
idx = find(alpha > Z);

% And update these chains
newgen(idx,1:MCMCPar.n+2) = [x(idx,1:MCMCPar.n) p_x(idx,1) log_p_x(idx,1)];

% And indicate that these chains have been accepted
accept(idx,1) = 1;