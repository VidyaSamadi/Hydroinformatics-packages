function [p,log_p,fx] = CompDensity(x,MCMCPar,Measurement,ModelName,Extra);
% This function computes the density of each x value

% Check whether to store the output of each model evaluation (function call)
if Measurement.N > 0,
    if strcmp(MCMCPar.modout,'Yes'),
        % Create initial ModPred
        fx = NaN(Measurement.N,size(x,1));
    end;
end;

% Loop over the individual parameter combinations of x
for ii = 1:size(x,1),

    % Call model to generate simulated data
    evalstr = ['fx(:,ii) = ',ModelName,'(x(ii,:),Extra);']; eval(evalstr);

    % If we have measured data --> calculate the residual (used by most likelihood functions)
    if Measurement.N > 0,
        % Calculate the error residual
        Err = (Measurement.MeasData(:) - fx(1:Measurement.N,ii));
    end;

    % Check whether the measurement error is estimated jointly with the parameters
    if isfield(Measurement,'Sigma') && isfield(Measurement,'n'),
        % If yes --> create the current value of Sigma from the inline function and ii-th candidate point 
        eval(Measurement.str_sigma);
    else
        % Sigma is either not needed in likelihood function or already defined numerically by the user
    end;

    % Model output is equal to posterior density (standard statistical distributions)
    if MCMCPar.lik == 1,
        p(ii,1) = fx(1,ii); log_p(ii,1) = log(p(ii,1));
    end;

    % Standard density function (iid error residuals)
    if MCMCPar.lik == 2,
        % Derive the log density
        if size(Sigma,1) == 1,  % --> homoscedastic error
            log_p(ii,1) = - ( Measurement.N / 2) * log(2 * pi) - Measurement.N * log( Sigma ) - 1/2 * Sigma^(-2) * sum ( Err.^2 );
        else                                % --> heteroscedastic error
            log_p(ii,1) = - ( Measurement.N / 2) * log(2 * pi) - sum ( log( Sigma ) ) - 1/2 * sum ( ( Err./Sigma ).^2);
        end;
        % And retain in memory
        p(ii,1) = log_p(ii,1);
    end;
        
    % same as previous one but now with sigma integrated out (derivation in lecture notes)
    if MCMCPar.lik == 3,
        % Derive the sum of squared error
        SSR = sum(abs(Err).^2);
        % And retain in memory
        p(ii,1) = -SSR; log_p(ii,1) = - Measurement.N/2 * log(SSR);
    end;

    % Model output is equal to log-posterior density
    if MCMCPar.lik == 4,
        p(ii,1) = fx(1,ii); log_p(ii,1) = p(ii,1);
    end;

    % Standard density function (iid error residuals) but with AR-1 model
    if MCMCPar.lik == 7,
        % First order autoregressive (AR-1) correction of residuals
        rho = x(ii,MCMCPar.n); Err_2 = Err(2:Measurement.N,1) - rho * Err(1:Measurement.N-1,1);
        % Now compute the log-likelihood
        if size(Sigma,1) == 1,  % --> homoscedastic error
            log_p(ii,1) = -(Measurement.N/2) * log(2*pi) - (Measurement.N/2) * log( Sigma^2 / (1-rho^2)) - ...
                (1/2) * (1-rho^2) * ( Err(1) / Sigma )^2 - (1/2) * sum( ( Err_2 ./ Sigma ).^2 );
        else                                % --> heteroscedastic error
            log_p(ii,1) = -(Measurement.N/2) * log(2*pi) - (Measurement.N/2) * log(mean( Sigma.^2 ) / (1-rho^2)) - ...
                (1/2) * (1-rho^2) * ( Err(1) / Sigma(1) )^2 - (1/2) * sum( ( Err_2 ./ Sigma(2:Measurement.N ) ).^2 );
        end;
        % And retain in memory
        p(ii,1) = log_p(ii,1);
    end;

    % Generalized log-likelihood (no Sigma needed)
    if MCMCPar.lik == 8,
        % Extract statistical model parameters
        par = Extra.fpar;               % fixed parameters
        par(Extra.idx_vpar) = x(ii,:);  % variable parameters
        par = par';                     % make it a column vector
        statpar = par(end-10:end);
        % Compute the log-likelihood
        log_p(ii,1) = GL('est',statpar,fx(1:Measurement.N,ii),Measurement.MeasData);
        % And retain in memory
        p(ii,1) = log_p(ii,1);
    end;

    % Whittle log-likelihood function
    if MCMCPar.lik == 9,
        % Calculate the log-likelihood using spectral density
        [log_L] = Whittle_logL(Measurement,fx(1:Measurement.N,ii));
        % Now store in memory
        p(ii,1) = log_L; log_p(ii,1) = p(ii,1);
    end;

    % Approximate Bayesian Computation (some estimate MCMCPar.delta along)
    if MCMCPar.lik == 10,
        % Now calculate rho
        rho = MCMCPar.rho( fx(1:Measurement.N,ii) , Measurement.MeasData(:) ) + normrnd(0,MCMCPar.delta,1,Measurement.N);
        % How many elements does rho consist of?
        N_rho = prod(size(rho));
        % Easier to work with log-density in practice --> when distance to 0 is large (with small delta)
        p(ii,1) = - ( N_rho / 2) * log(2 * pi) - N_rho * log( MCMCPar.delta ) - 1/2 * MCMCPar.delta^(-2) * sum ( rho.^2 );
        % And log density
        log_p(ii,1) = p(ii,1);
    end;

    % Approximate Bayesian Computation (alternative to continuous kernel)
    if MCMCPar.lik == 11,
        % Now calculate rho
        p(ii,1) = max ( abs ( MCMCPar.rho( fx(1:Measurement.N,ii) , Measurement.MeasData(:) ) ) );
        % And log density
        log_p(ii,1) = p(ii,1);
    end;

end; 