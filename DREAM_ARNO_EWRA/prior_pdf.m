function [log_prior] = prior_pdf(MCMCPar,x);
% Calculate the alpha ratio of the prior densities

% How many members does x have?
N = size(x,1);

% Check whether an explicit prior distribution is used
if strcmp(lower(MCMCPar.prior),'prior');

    % Initialize prior_x
    prior_x = NaN(N,MCMCPar.n);
    
    % Compute prior densities for each parameter in each chain
    for qq = 1 : MCMCPar.n,
        % And for each proposal
        for zz = 1 : N,
            % Compute prior of proposal
            prior_x(zz,qq) = eval(char(strrep(MCMCPar.prior_marginal(qq),'rnd(','pdf(x(zz,qq),')));
        end;
    end;

    % Calculate the prior density (product of individual parameter densities)
    log_prior = sum( log ( max ( prior_x , 1e-299 ) ) , 2);
    
else
    
    % Simply set to 0
    log_prior = zeros(N,1);
    
end;