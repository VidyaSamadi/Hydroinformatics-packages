clc
clear
% Function runs GLUE
Model = 'Swap';

if strcmp(Model,'Swap');

    MCMCPar.n = 5;                     % Dimension of the problem
    MCMCPar.ndraw = 10000;              % Maximum number of function evaluations
     
    % Define feasible parameter space (minimum and maximum values)
    %  parameters:      ORES, OSAT, ALFA, NPAR, KSAT, LEXP
        ParRange.minn = [1 0.1 0.1 0 0.1];
        ParRange.maxn = [500 2 0.99 0.1 0.99];
        
    % Latin Hypercube sampling (options, "LHS", and "PRIOR") 
        MCMCPar.prior = 'LHS';
 
    % Check Nsim in Glue.m at line 71,note that it should be equal
    % with numbers of observations
        
    % Store in array Measurement.MeasData
    load bound.txt;
    MaxT = size(bound,1);
    Measurement.MeasData = bound(1:MaxT,4);
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
save output