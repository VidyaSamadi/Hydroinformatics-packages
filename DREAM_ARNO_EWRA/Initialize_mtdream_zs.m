function [MCMCPar,chain,X,output,fx,fid_Z,delta_tot,pCR,lCR,...
    CR,Iter,iteration,T,T2,Table_JumpRate,iloc] = Initialize_mtdream_zs(MCMCPar,Func_name,Extra,ParRange,Measurement);
% Initializes the main variables used in DREAM_ZS

% Check whether to do a restart or not?
if strcmp(lower(MCMCPar.Restart),'no')
    
    % Set random number generator
    opts.Seed = 'sum(100*clock)  % evaluated if it is a string';
    % Then generate new seed
    if ischar(opts.Seed)
        randn('state', eval(opts.Seed));     % random number generator state
    else
        randn('state', opts.Seed);
    end
    
    % Initialize basic variables
    [MCMCPar,pCR,lCR,CR,Iter,iteration,T,T2,chain,Table_JumpRate,iloc,output] = InitMCMC(MCMCPar,Measurement);
    
    % Step 1: Sample MCMCPar.m0 points in the parameter space and store in Z
    if strcmp(lower(MCMCPar.prior),'lhs'),
        % Initialize chains with Latin hypercube sampling
        [Zinit] = LHS(ParRange.minn,ParRange.maxn,MCMCPar.m0 + MCMCPar.seq);
    elseif strcmp(lower(MCMCPar.prior),'cov');
        % Initialize chains with (multi)-normal distribution
        [Zinit] = repmat(MCMCPar.mu,MCMCPar.m0 + MCMCPar.seq,1) + randn(MCMCPar.m0 + MCMCPar.seq,MCMCPar.n) * chol(MCMCPar.cov);
    elseif strcmp(lower(MCMCPar.prior),'grid');
        % Initialize chains with grid distribution
        F = [0 : 1/(MCMCPar.m0 + MCMCPar.seq - 1) : 1]'; ii = randperm ( MCMCPar.m0 + MCMCPar.seq ); F = F(ii);
        % Now create the grid (useful to avoid problems with density functions ==> 0)
        for qq = 1 : MCMCPar.m0 + MCMCPar.seq,
            % Create vectors sequentially
            [Zinit(qq,1:MCMCPar.n)] = ParRange.minn + F(qq) * ( ParRange.maxn - ParRange.minn );
        end;
    elseif strcmp(lower(MCMCPar.prior),'prior');
        % Create the initial position of each chain by drawing each parameter individually from the prior
        for qq = 1:MCMCPar.n,
            for zz = 1:MCMCPar.m0 + MCMCPar.seq,
                Zinit(zz,qq) = eval(char(MCMCPar.prior_marginal(qq)));
            end;
        end;
    end;
    
    % Do boundary handling -- what to do when points fall outside bound
    if strcmp(lower(MCMCPar.BoundHandling),'none') == 0;
        [Zinit] = BoundaryHandling(Zinit,ParRange,MCMCPar.BoundHandling);
    end;
    
    % Define initial MCMCPar.m0 rows of Z to be initial sample -- posterior density is not needed and thus not evaluated!!
    Z(1:MCMCPar.m0,1:MCMCPar.n+2) = [ Zinit(1:MCMCPar.m0,1:MCMCPar.n) NaN(MCMCPar.m0,2) ];
    
    % Now write Z to binary file
    fid_Z = fopen('Z.bin','w+','n');
    
    % Now write --> Z' is required to be able to use fseek later on
    fwrite(fid_Z,Z','double');
    
    % Now close file
    fclose(fid_Z);
    
    % Now open the file to read
    fid_Z = fopen('Z.bin','r');
    
    % Define initial population from last MCMCPar.seq samples of Zinit
    X = Zinit(MCMCPar.m0 + 1:MCMCPar.m0+MCMCPar.seq,1:MCMCPar.n); clear Zinit;
    
    % Now evaluate the model ( = pdf ) and return fx
    [fx(:,1:MCMCPar.seq)] = Evaluate_model(X,MCMCPar,Measurement,Func_name,Extra);
    
    % Calculate the likelihood and log-likelihood from fx
    [p,log_p] = CalcLikelihood(X,fx(:,1:MCMCPar.seq),MCMCPar,Measurement);
    
    % Append X with information about posterior density (or transformation thereof)
    X = [X p log_p];
    
    % Store the model simulations (if appropriate)
    [ dummy ] = mtdream_zs_store_results ( MCMCPar , fx , Measurement , 'w+' );

    % Set the first point of each of the MCMCPar.seq chain equal to the initial X values
    chain(1,1:MCMCPar.n+2,1:MCMCPar.seq) = reshape(X',1,MCMCPar.n+2,MCMCPar.seq);
    
    % Save N_CR in memory and initialize sum_p2
    output.CR(1,1:MCMCPar.nCR+1) = [Iter pCR]; delta_tot = zeros(1,MCMCPar.nCR);
    
    % Compute the R-statistic of Gelman and Rubin
    [output.R_stat(1,1:MCMCPar.n+1)] = [Iter Gelman(chain(1:iloc,1:MCMCPar.n,1:MCMCPar.seq),MCMCPar)];
    
elseif strcmp(lower(MCMCPar.Restart),'yes'),
    
    % If a restart run is being done: just load the output from the previous ongoing trial
    load MT-DREAM_ZS.mat;

    % Open external archive as well
    fid_Z = fopen('Z.bin','r');

    % And make sure we add zeros to "chain" array
    chain = [chain ; zeros(size(chain,1)-1,size(chain,2),size(chain,3))];
    
    % And double population
    MCMCPar.ndraw = 2 * MCMCPar.ndraw;
    
end;