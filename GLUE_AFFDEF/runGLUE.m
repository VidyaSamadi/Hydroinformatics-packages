% Function runs GLUE
Model = 'Swap';

if strcmp(Model,'Swap');

    MCMCPar.n = 6;                     % Dimension of the problem
    MCMCPar.ndraw = 50;              % Maximum number of function evaluations
     
    % Define feasible parameter space (minimum and maximum values)
    %  parameters:      ORES, OSAT, ALFA, NPAR, KSAT, LEXP
        ParRange.minn = [0.04 0.42 0.008 1.10 5.0 -6.1];
        ParRange.maxn = [0.15 0.53 0.040 1.57 50.0 -1.1];
        
    % Latin Hypercube sampling (options, "LHS", and "PRIOR") 
        MCMCPar.prior = 'PRIOR';
 
    % Check Nsim in Glue.m at line 71,note that it should be equal
    % with numbers of observations
        
    % Store in array Measurement.MeasData
    load obs_moist.txt;
    MaxT = size(obs_moist,1);
    Measurement.MeasData = obs_moist(1:MaxT,1);
    Measurement.N = size(Measurement.MeasData,1);
    
    % Specify the prior distributions for the various parameters
    MCMCPar.prior_marginal = {  'normrnd(0.0973,0.0128)',...
                                'normrnd(0.4760,0.0137)',...
                                'normrnd(-1.818,0.108)',...
                                'normrnd(0.1220,0.0189)',...
                                'normrnd(1.442,0.161)',...
                                'normrnd(-1.1760,1.196)'};
                              
    
    % Define modelName
    ModelName = 'Swap';
    
    % Run GLUE
    [output] = GLUE(MCMCPar,ModelName,Measurement,ParRange);
    
end;