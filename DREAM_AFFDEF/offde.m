function [x_new,CR] = offde(x_old,X,CR,MCMCPar,Table_JumpRate,ParRange,BoundHandling,R2,DR)
%Generates offspring using METROPOLIS HASTINGS monte-carlo markov chain

% Generate ergodicity term
eps = 1e-6 * randn(MCMCPar.seq,MCMCPar.n);

% If not a delayed rejection step --> generate proposal with DE
if strcmp(DR,'No');

    % Determine which sequences to evolve with what DE strategy
    [DEversion] = DEStrategy(MCMCPar);

    % Generate series of permutations of chains
    [dummy,tt] = sort(rand(MCMCPar.seq-1,MCMCPar.seq));

    % Generate uniform random numbers for each chain to determine which dimension to update
    D = rand(MCMCPar.seq,MCMCPar.n);

    % Ergodicity for each individual chain
    noise_x = MCMCPar.eps * (2 * rand(MCMCPar.seq,MCMCPar.n) - 1);

    % Initialize the delta update to zero
    delta_x = zeros(MCMCPar.seq,MCMCPar.n);

    % Each chain evolves using information from other chains to create offspring
    for qq = 1:MCMCPar.seq,

        % Define ii and remove current member as an option
        ii = ones(MCMCPar.seq,1); ii(qq) = 0; idx = find(ii>0);

        % randomly select two members of ii that have value == 1
        rr = idx(tt(1:2*DEversion(qq,1),qq));

        % --- WHICH DIMENSIONS TO UPDATE? DO SOMETHING WITH CROSSOVER ----
        [i] = find(D(qq,1:MCMCPar.n) > (1-CR(qq,1)));

        % Update at least one dimension
        if isempty(i), i = randperm(MCMCPar.n); i = i(1); end;
        % ----------------------------------------------------------------

        % Determine the number of dimensions that are going to be updated
        NrDim = size(i,2);
        
        % Determine the associated JumpRate and compute the jump
        if (rand < 4/5),
            % Lookup Table
            JumpRate = Table_JumpRate(NrDim,DEversion(qq,1));
            
            % Produce the difference of the pairs used for population evolution
            delta = sum(X(rr(1:DEversion(qq,1)),1:MCMCPar.n) - X(rr(DEversion(qq,1)+1:2*DEversion(qq,1)),1:MCMCPar.n),1);
            
            % Then fill update the dimension
            delta_x(qq,i) = (1 + noise_x(qq,i)) * JumpRate.*delta(1,i);
        else
            
            % Set the JumpRate to 1 and overwrite CR and DEversion
            JumpRate = 1; CR(qq,1) = -1; 
            
            % Compute delta from one pair
            delta = X(rr(1),1:MCMCPar.n) - X(rr(2),1:MCMCPar.n);
            
            % Now jumprate to facilitate jumping from one mode to the other in all dimensions
            delta_x(qq,1:MCMCPar.n) = JumpRate * delta;
        end;
        
        % Check this line to avoid that jump = 0 and xnew is similar to xold
        if (sum(delta_x(qq,1:MCMCPar.n).^2,2) == 0),
            % Compute the Cholesky Decomposition of x_old
            R = (2.38/sqrt(MCMCPar.n)) * chol(cov(x_old(1:end,1:MCMCPar.n)) + 1e-5*eye(MCMCPar.n));           
            % Generate jump using multinormal distribution
            delta_x(qq,1:MCMCPar.n) = randn(1,MCMCPar.n) * R;
        end;
        
    end;
end;

% If delayed rejection step --> generate proposal with DR
if strcmp(DR,'Yes');
    % Loop over all chains -- all dimensions are updated
    for qq = 1:MCMCPar.seq,
        % Generate a new proposal distance using standard procedure
        delta_x(qq,1:MCMCPar.n) = randn(1,MCMCPar.n) * R2;
    end;
end;

% Update x_old with delta_x and eps;
x_new = x_old + delta_x + eps;

% Do boundary handling -- what to do when points fall outside bound
if strcmp(BoundHandling,'Reflect');
    [x_new] = ReflectBounds(x_new,ParRange);
end;
if strcmp(BoundHandling,'Bound');
    [x_new] = SetToBounds(x_new,ParRange);
end;
if strcmp(BoundHandling,'Fold');
    [x_new] = FoldBounds(x_new,ParRange);
end;
if strcmp(BoundHandling,'None');
    % Do nothing
end;